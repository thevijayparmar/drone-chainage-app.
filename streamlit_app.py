import streamlit as st
import re, math, zipfile, io, xml.etree.ElementTree as ET
import plotly.graph_objects as go

st.set_page_config(page_title="Drone SRT & Chainage", layout="wide")
st.title("Drone SRT & Chainage Workflow")

# -- Helpers --------------------------------------------------------------

def hav(p, q):
    R = 6371.0088
    φ1, φ2 = math.radians(p[0]), math.radians(q[0])
    dφ = math.radians(q[0] - p[0]); dλ = math.radians(q[1] - p[1])
    a = math.sin(dφ/2)**2 + math.cos(φ1)*math.cos(φ2)*math.sin(dλ/2)**2
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))

@st.cache_data
def parse_srt(uploaded):
    lines = uploaded.getvalue().decode('utf-8').splitlines()
    blocks=[]; i=0
    while i < len(lines):
        if not lines[i].isdigit(): i+=1; continue
        idx=int(lines[i]); rng=lines[i+1]; j=i+2
        lat=lon=None; alt=0.0; tim=""
        while j < len(lines) and lines[j].strip():
            L=lines[j]
            m=re.search(r"\[latitude:\s*([0-9.\-]+)",L)
            if m: lat=float(m.group(1))
            m=re.search(r"\[longitude:\s*([0-9.\-]+)",L)
            if m: lon=float(m.group(1))
            m=re.search(r"\[altitude:\s*([0-9.\-]+)",L)
            if m: alt=float(m.group(1))
            m=re.search(r"\d{4}-\d{2}-\d{2}\s*([0-9:]{8})",L)
            if m: tim=m.group(1)
            j+=1
        if lat is not None and lon is not None:
            blocks.append({ 'idx':idx, 'range':rng, 'lat':lat, 'lon':lon, 'alt':alt, 'tim':tim })
        i=j+1
    return blocks

@st.cache_data
def parse_kml_2d(uploaded):
    tree = ET.parse(io.BytesIO(uploaded.getvalue()))
    ns={'kml':'http://www.opengis.net/kml/2.2'}
    coords=[]
    for e in tree.findall('.//kml:LineString/kml:coordinates', ns):
        for c in e.text.strip().split():
            lon,lat,*_ = map(float,c.split(','))
            coords.append((lat, lon))
    return coords

@st.cache_data
def cumulative_dist(coords):
    cum=[0.0]
    for A,B in zip(coords, coords[1:]): cum.append(cum[-1]+hav(A,B))
    return cum

# markers every 50m
@st.cache_data
def generate_markers(coords, cumd, interval_km=0.05):
    markers=[]; d=0.0; idx=0; total=cumd[-1]
    while d<=total:
        while idx<len(cumd)-1 and cumd[idx+1]<d: idx+=1
        A,B=coords[idx],coords[idx+1]
        seg=cumd[idx+1]-cumd[idx]
        f=(d-cumd[idx])/seg if seg>0 else 0
        markers.append((A[0]+f*(B[0]-A[0]), A[1]+f*(B[1]-A[1]), d))
        d+=interval_km
    return markers

# build chainage KML string
@st.cache_data
def kml_from_markers(markers):
    lines=["<?xml version='1.0' encoding='UTF-8'?>",
           "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>"]
    for lat,lon,d in markers:
        lines.append(f"<Placemark><name>{d:.3f} km</name><Point><coordinates>{lon},{lat},0</coordinates></Point></Placemark>")
    lines.append("</Document></kml>")
    return '\n'.join(lines)

# project onto chain
@st.cache_data
def project_to_line(pt, coords, cumd):
    dists=[hav(pt,c) for c in coords]
    i=min(range(len(dists)), key=lambda i:dists[i])
    return i, cumd[i]

# build SRT route KML
@st.cache_data
def kml_from_srt(blocks, cumd, si, ei):
    lines=["<?xml version='1.0' encoding='UTF-8'?>",
           "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>",
           "<Placemark><LineString><coordinates>"]
    for b in blocks: lines.append(f"{b['lon']},{b['lat']},{b['alt']}")
    lines.append("</coordinates></LineString></Placemark>")
    s,e=blocks[si],blocks[ei]
    lines.append(f"<Placemark><name>START {cumd[si]:.3f} km</name><Point><coordinates>{s['lon']},{s['lat']},{s['alt']}</coordinates></Point></Placemark>")
    lines.append(f"<Placemark><name>END {cumd[ei]:.3f} km</name><Point><coordinates>{e['lon']},{e['lat']},{e['alt']}</coordinates></Point></Placemark>")
    lines.append("</Document></kml>")
    return '\n'.join(lines)

# write SRTs to zip
@st.cache_data
def srt_zip(blocks, prefix):
    buf=io.BytesIO(); z=zipfile.ZipFile(buf,'w')
    cum=cumulative_dist([(b['lat'],b['lon']) for b in blocks])
    mapping={
        f"{prefix}_01_Latitude.srt":[f"{b['lat']:.6f}" for b in blocks],
        f"{prefix}_02_Longitude.srt":[f"{b['lon']:.6f}" for b in blocks],
        f"{prefix}_03_Altitude.srt":[f"{b['alt']:.2f}" for b in blocks],
        f"{prefix}_04_Timer.srt":[b['tim'] for b in blocks],
        f"{prefix}_05_LatLonOutput.srt":[f"{b['lat']:.6f}, {b['lon']:.6f}" for b in blocks],
        f"{prefix}_06_Chainage.srt":[f"{d:.3f} km" for d in cum]
    }
    for fn,txts in mapping.items():
        content='\n\n'.join(f"{blocks[i]['idx']}\n{blocks[i]['range']}\n{txts[i]}" for i in range(len(blocks)))
        z.writestr(fn, content)
    z.close(); buf.seek(0)
    return buf.getvalue()

# -- UI ------------------------------------------------------------------

# Step 1: upload
col1, col2 = st.columns(2)
with col1: up_srt = st.file_uploader("DJI .srt file", type='srt')
with col2: up_kml = st.file_uploader("Full-chainage KML", type=['kml','xml'])
if up_srt and up_kml:
    chain_coords = parse_kml_2d(up_kml)
    cum_chain    = cumulative_dist(chain_coords)

    # Step 2+3: generate markers & confirm
    if 'chain_rev' not in st.session_state: st.session_state.chain_rev=False
    if st.button("▶ Generate & Download 50 m Markers KML"):
        markers = generate_markers(chain_coords[::-1] if st.session_state.chain_rev else chain_coords, cum_chain[::-1] if st.session_state.chain_rev else cum_chain)
        kml = kml_from_markers(markers)
        st.download_button("Download markers_50m.kml", kml, file_name="markers_50m.kml", mime="application/vnd.google-earth.kml+xml")
    dir_ok = st.radio("Are 50 m markers direction correct?", ["Yes","Reverse Needed"])
    if dir_ok=="Reverse Needed":
        st.session_state.chain_rev = not st.session_state.chain_rev

    # proceed only if confirmed
    if dir_ok=="Yes":
        if st.session_state.chain_rev:
            chain_coords.reverse(); cum_chain.reverse()

        # Steps 4–5: parse SRT & project
        blocks0 = parse_srt(up_srt)
        coords0 = [(b['lat'],b['lon']) for b in blocks0]
        alts0   = [b['alt'] for b in blocks0]
        cum0    = cumulative_dist(coords0)
        si0,km_s0 = project_to_line(coords0[0], chain_coords, cum_chain)
        ei0,km_e0 = project_to_line(coords0[-1],chain_coords, cum_chain)
        st.write(f"Start → {km_s0:.3f} km, End → {km_e0:.3f} km, Total → {cum0[-1]*1000:.1f} m")
        dir2 = st.radio("Is SRT start→end correct?", ["Yes","Reverse Needed"], key="srtdir")
        if dir2=="Reverse Needed":
            blocks0.reverse(); coords0.reverse(); alts0.reverse(); cum0 = cumulative_dist(coords0)
            si0,km_s0 = project_to_line(coords0[0], chain_coords, cum_chain)
            ei0,km_e0 = project_to_line(coords0[-1],chain_coords, cum_chain)

        # Step 6: 3D plot
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x=[c[1] for c in chain_coords], y=[c[0] for c in chain_coords], z=[0]*len(chain_coords), line=dict(color='gray'), name='Chainage'))
        fig.add_trace(go.Scatter3d(x=[c[1] for c in coords0], y=[c[0] for c in coords0], z=alts0, line=dict(color='red'), name='SRT'))
        fig.add_trace(go.Scatter3d(x=[coords0[0][1]], y=[coords0[0][0]], z=[alts0[0]], mode='markers+text', text=[f"Start {km_s0:.3f} km"], marker=dict(size=5,color='green'), name='Start'))
        fig.add_trace(go.Scatter3d(x=[coords0[-1][1]], y=[coords0[-1][0]], z=[alts0[-1]], mode='markers+text', text=[f"End {km_e0:.3f} km"], marker=dict(size=5,color='blue'), name='End'))
        fig.update_layout(scene=dict(xaxis_title='Lon', yaxis_title='Lat', zaxis_title='Alt (m)'), height=600)
        st.plotly_chart(fig)

        # Step 7: export initial
        kml0 = kml_from_srt(blocks0, cum0, si0, ei0)
        st.download_button("Download initial_srt_route.kml", kml0, "initial_srt_route.kml", mime="application/vnd.google-earth.kml+xml")
        zip0 = srt_zip(blocks0, "initial")
        st.download_button("Download initial SRT files ZIP", zip0, "initial_srt_files.zip", mime="application/zip")

        # Step 8: bulk
        up_bulk = st.file_uploader("Upload additional .srt files (bulk)", type='srt', accept_multiple_files=True)
        if up_bulk:
            all_buf = io.BytesIO(); zall=zipfile.ZipFile(all_buf,'w')
            for bf in up_bulk:
                name = bf.name[:-4]
                blks = parse_srt(bf)
                crd = [(b['lat'],b['lon']) for b in blks]
                altb= [b['alt'] for b in blks]
                cumb= cumulative_dist(crd)
                si,ks = project_to_line(crd[0], chain_coords, cum_chain)
                ei,ke = project_to_line(crd[-1],chain_coords, cum_chain)
                kmlb = kml_from_srt(blks, cumb, si, ei)
                zall.writestr(f"{name}_route.kml", kmlb)
                for fn in write_srt_files_with_prefix(blks,name):
                    zall.write(fn)
            zall.close(); all_buf.seek(0)
            st.download_button("Download all bulk outputs ZIP", all_buf.read(), "all_bulk_outputs.zip", mime="application/zip")

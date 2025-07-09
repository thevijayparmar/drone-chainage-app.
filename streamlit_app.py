import streamlit as st
import re, math, zipfile, io, xml.etree.ElementTree as ET
import plotly.graph_objects as go

# ──────────────────────────────────────────────────────────────────────────
# Author: Vijay Parmar
# Community: BGol Community of Advanced Surveying and GIS Professionals
# ──────────────────────────────────────────────────────────────────────────

st.set_page_config(page_title="Drone SRT & Chainage", layout="wide")
st.title("Drone SRT & Chainage Workflow")
st.markdown(
    """**Author:** Vijay Parmar  
**Community:** BGol Community of Advanced Surveying and GIS Professionals"""
)

# ──────────────────────────────────────────────────────────────────────────
# Helper Functions
# ──────────────────────────────────────────────────────────────────────────

def hav(p, q):
    R = 6371.0088
    φ1, φ2 = math.radians(p[0]), math.radians(q[0])
    dφ = math.radians(q[0] - p[0]); dλ = math.radians(q[1] - p[1])
    a = math.sin(dφ/2)**2 + math.cos(φ1)*math.cos(φ2)*math.sin(dλ/2)**2
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))

@st.cache_data
def parse_srt(uploaded):
    lines = uploaded.getvalue().decode('utf-8').splitlines()
    blocks, i = [], 0
    while i < len(lines):
        if not lines[i].strip().isdigit():
            i += 1
            continue
        idx = int(lines[i].strip())
        rng = lines[i+1].strip()
        j = i + 2
        lat = lon = None
        alt = 0.0
        tim = ""
        while j < len(lines) and lines[j].strip():
            L = lines[j]
            if m := re.search(r"\[latitude:\s*([0-9.\-]+)", L):
                lat = float(m.group(1))
            if m := re.search(r"\[longitude:\s*([0-9.\-]+)", L):
                lon = float(m.group(1))
            if m := re.search(r"\[altitude:\s*([0-9.\-]+)", L):
                alt = float(m.group(1))
            if m := re.search(r"\d{4}-\d{2}-\d{2}\s*([0-9:]{8})", L):
                tim = m.group(1)
            j += 1
        if lat is not None and lon is not None:
            blocks.append({
                'idx': idx, 'range': rng,
                'lat': lat, 'lon': lon,
                'alt': alt, 'tim': tim
            })
        i = j + 1
    return blocks

@st.cache_data
def parse_kml_2d(uploaded):
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    tree = ET.parse(io.BytesIO(uploaded.getvalue()))
    coords = []
    for e in tree.findall('.//kml:LineString/kml:coordinates', ns):
        for part in e.text.strip().split():
            lon, lat, *_ = part.split(',')
            coords.append((float(lat), float(lon)))
    return coords

@st.cache_data
def cumulative_dist(coords):
    cum = [0.0]
    for A, B in zip(coords, coords[1:]):
        cum.append(cum[-1] + hav(A, B))
    return cum

@st.cache_data
def generate_markers(coords, cumd, interval_km=0.05):
    markers, d, idx = [], 0.0, 0
    total = cumd[-1]
    while d <= total:
        while idx < len(cumd) - 1 and cumd[idx+1] < d:
            idx += 1
        A, B = coords[idx], coords[idx+1]
        seg = cumd[idx+1] - cumd[idx]
        frac = (d - cumd[idx]) / seg if seg > 0 else 0
        markers.append((
            A[0] + frac*(B[0] - A[0]),
            A[1] + frac*(B[1] - A[1]),
            d
        ))
        d += interval_km
    return markers

@st.cache_data
def kml_from_markers(markers):
    lines = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        "<!-- Author: Vijay Parmar; Community: BGol Community of Advanced Surveying and GIS Professionals -->",
        "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>"
    ]
    for lat, lon, d in markers:
        lines.append(
            f"<Placemark>"
            f"<name>{d:.3f} km</name>"
            f"<Point><coordinates>{lon:.6f},{lat:.6f},0</coordinates></Point>"
            f"</Placemark>"
        )
    lines.append("</Document></kml>")
    return "\n".join(lines)

@st.cache_data
def project_to_line(pt, coords, cumd):
    dists = [hav(pt, c) for c in coords]
    i = min(range(len(dists)), key=lambda i: dists[i])
    return i, cumd[i]

@st.cache_data
def kml_from_srt(blocks, cumd, si, ei):
    lines = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        "<!-- Author: Vijay Parmar; Community: BGol Community of Advanced Surveying and GIS Professionals -->",
        "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>",
        "<Style id='routeStyle'>",
        "  <LineStyle><color>ff0000ff</color><width>10.0</width></LineStyle>",
        "</Style>"
    ]
    lines.append("<Placemark><styleUrl>#routeStyle</styleUrl><name>SRT Route</name><LineString><coordinates>")
    for b in blocks:
        lines.append(f"{b['lon']:.6f},{b['lat']:.6f},{b['alt']:.2f}")
    lines.append("</coordinates></LineString></Placemark>")
    s, e = blocks[si], blocks[ei]
    lines.append(
        f"<Placemark><styleUrl>#routeStyle</styleUrl>"
        f"<name>START {cumd[si]:.3f} km</name>"
        f"<Point><coordinates>{s['lon']:.6f},{s['lat']:.6f},{s['alt']:.2f}</coordinates></Point>"
        f"</Placemark>"
    )
    lines.append(
        f"<Placemark><styleUrl>#routeStyle</styleUrl>"
        f"<name>END {cumd[ei]:.3f} km</name>"
        f"<Point><coordinates>{e['lon']:.6f},{e['lat']:.6f},{e['alt']:.2f}</coordinates></Point>"
        f"</Placemark>"
    )
    lines.append("</Document></kml>")
    return "\n".join(lines)

@st.cache_data
def srt_zip(blocks, prefix, offset=0.0):
    buf = io.BytesIO()
    z = zipfile.ZipFile(buf, 'w')
    cum = cumulative_dist([(b['lat'], b['lon']) for b in blocks])
    cum = [offset + d for d in cum]
    mapping = {
        f"{prefix}_01_Latitude.srt":     [f"{b['lat']:.6f}" for b in blocks],
        f"{prefix}_02_Longitude.srt":    [f"{b['lon']:.6f}" for b in blocks],
        f"{prefix}_03_Altitude.srt":     [f"{b['alt']:.2f}" for b in blocks],
        f"{prefix}_04_Timer.srt":        [b['tim'] for b in blocks],
        f"{prefix}_05_LatLonOutput.srt": [f"{b['lat']:.6f}, {b['lon']:.6f}" for b in blocks],
        f"{prefix}_06_Chainage.srt":     [f"{d:.3f} km" for d in cum]
    }
    for fn, txts in mapping.items():
        content = "\n\n".join(
            f"{blocks[i]['idx']}\n{blocks[i]['range']}\n{txts[i]}"
            for i in range(len(blocks))
        )
        z.writestr(fn, content)
    z.close(); buf.seek(0)
    return buf.read()

# ──────────────────────────────────────────────────────────────────────────
# 1. Upload Inputs
# ──────────────────────────────────────────────────────────────────────────

col1, col2 = st.columns(2)
with col1:
    up_srt = st.file_uploader("DJI .srt file", type='srt')
with col2:
    up_kml = st.file_uploader("Full-chainage .kml/.xml", type=['kml','xml'])

if up_srt and up_kml:
    base_chain = parse_kml_2d(up_kml)

    # Step 2: markers + FLIP + proceed
    st.subheader("Step 2: Generate 50 m Markers")
    flip = st.checkbox("FLIP chainage direction", key="flip_chain")
    chain_coords = list(reversed(base_chain)) if flip else base_chain
    cum_chain    = cumulative_dist(chain_coords)

    if st.button("▶ Generate & Download 50 m Markers KML"):
        kmlm = kml_from_markers(generate_markers(chain_coords, cum_chain))
        st.download_button(
            "Download markers_50m.kml", kmlm, "markers_50m.kml",
            mime="application/vnd.google-earth.kml+xml"
        )

    if "proceed" not in st.session_state:
        st.session_state.proceed = False
    if st.button("Proceed further"):
        st.session_state.proceed = True

    # Steps 4–9
    if st.session_state.proceed:
        # Step 4: project start/end
        blocks0 = parse_srt(up_srt)
        coords0 = [(b['lat'],b['lon']) for b in blocks0]
        alts0   = [b['alt'] for b in blocks0]
        cum0    = cumulative_dist(coords0)
        si0, sk0 = project_to_line(coords0[0], chain_coords, cum_chain)
        ei0, ek0 = project_to_line(coords0[-1], chain_coords, cum_chain)

        st.subheader("Step 4: SRT Stationing")
        st.write(f"• Start → **{sk0:.3f} km**")
        st.write(f"• End   → **{ek0:.3f} km**")
        st.write(f"• Total flight length: **{cum0[-1]*1000:.1f} m**")

        # Step 6: 3D visualize
        st.subheader("Step 6: 3D Visualization")
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(
            x=[c[1] for c in chain_coords],
            y=[c[0] for c in chain_coords],
            z=[0]*len(chain_coords),
            line=dict(color='blue', width=10.0),
            name='Chainage'
        ))
        fig.add_trace(go.Scatter3d(
            x=[c[1] for c in coords0],
            y=[c[0] for c in coords0],
            z=alts0,
            line=dict(color='red', width=4.0),
            name='SRT Route'
        ))
        fig.add_trace(go.Scatter3d(
            x=[coords0[0][1]], y=[coords0[0][0]], z=[alts0[0]],
            mode='markers+text',
            text=[f"Start {sk0:.3f} km"],
            marker=dict(size=5,color='green'),
            name='Start'
        ))
        fig.add_trace(go.Scatter3d(
            x=[coords0[-1][1]], y=[coords0[-1][0]], z=[alts0[-1]],
            mode='markers+text',
            text=[f"End {ek0:.3f} km"],
            marker=dict(size=5,color='blue'),
            name='End'
        ))
        fig.update_layout(scene=dict(
            xaxis_title='Lon', yaxis_title='Lat', zaxis_title='Alt (m)'
        ), height=600)
        st.plotly_chart(fig)

        # Step 7: export initial
        st.subheader("Step 7: Export Initial SRT‐Derived Files")
        kml0 = kml_from_srt(blocks0, [sk0 + d for d in cum0], si0, ei0)
        st.download_button("Download initial_srt_route.kml", kml0,
                           "initial_srt_route.kml",
                           mime="application/vnd.google-earth.kml+xml")
        zip0 = srt_zip(blocks0, "initial", offset=sk0)
        st.download_button("Download initial SRT files ZIP", zip0,
                           "initial_srt_files.zip", mime="application/zip")

        # Step 8: bulk
        st.subheader("Step 8: Bulk Processing")
        up_bulk = st.file_uploader(
            "Upload additional .srt files (bulk)",
            type='srt', accept_multiple_files=True
        )
        if up_bulk:
            all_buf = io.BytesIO()
            zall = zipfile.ZipFile(all_buf, 'w')
            for bf in up_bulk:
                name = bf.name.rsplit('.',1)[0]
                blks = parse_srt(bf)
                crd  = [(b['lat'], b['lon']) for b in blks]
                altb = [b['alt'] for b in blks]
                cumb = cumulative_dist(crd)
                si, sk = project_to_line(crd[0], chain_coords, cum_chain)
                ei, ek = project_to_line(crd[-1], chain_coords, cum_chain)
                kmlb = kml_from_srt(blks, [sk + d for d in cumb], si, ei)
                zall.writestr(f"{name}_route.kml", kmlb)
                for fn in srt_zip(blks, name, offset=sk):
                    zall.write(fn)
            zall.close(); all_buf.seek(0)
            st.download_button("Download all bulk outputs ZIP", all_buf.read(),
                               "all_bulk_outputs.zip", mime="application/zip")

        st.success("✅ All outputs are ready for download.")

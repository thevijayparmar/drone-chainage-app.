import streamlit as st
import re, math, zipfile, io, xml.etree.ElementTree as ET
import plotly.graph_objects as go

st.set_page_config(page_title="Drone SRT & Chainage", layout="wide")
st.title("Drone SRT & Chainage Workflow")

# -- Helper Functions -----------------------------------------------------

def hav(p, q):
    R = 6371.0088
    φ1, φ2 = math.radians(p[0]), math.radians(q[0])
    dφ = math.radians(q[0] - p[0])
    dλ = math.radians(q[1] - p[1])
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
            if m := re.search(r"\[latitude:\s*([0-9.\-]+)", L): lat = float(m.group(1))
            if m := re.search(r"\[longitude:\s*([0-9.\-]+)", L): lon = float(m.group(1))
            if m := re.search(r"\[altitude:\s*([0-9.\-]+)", L): alt = float(m.group(1))
            if m := re.search(r"\d{4}-\d{2}-\d{2}\s*([0-9:]{8})", L): tim = m.group(1)
            j += 1
        if lat is not None and lon is not None:
            blocks.append({ 'idx': idx, 'range': rng, 'lat': lat, 'lon': lon, 'alt': alt, 'tim': tim })
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
    markers = []
    total = cumd[-1]
    d = 0.0
    idx = 0
    while d <= total:
        while idx < len(cumd) - 1 and cumd[idx+1] < d:
            idx += 1
        A, B = coords[idx], coords[idx+1]
        seg = cumd[idx+1] - cumd[idx]
        frac = (d - cumd[idx]) / seg if seg > 0 else 0
        lat = A[0] + frac * (B[0] - A[0])
        lon = A[1] + frac * (B[1] - A[1])
        markers.append((lat, lon, d))
        d += interval_km
    return markers

@st.cache_data
def kml_from_markers(markers):
    lines = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>"
    ]
    for lat, lon, d in markers:
        lines.append(
            f"<Placemark><name>{d:.3f} km</name>"
            f"<Point><coordinates>{lon},{lat},0</coordinates></Point></Placemark>"
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
        "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>",
        "<Placemark><name>SRT Route</name><LineString><coordinates>"
    ]
    for b in blocks:
        lines.append(f"{b['lon']:.6f},{b['lat']:.6f},{b['alt']:.2f}")
    lines.append("</coordinates></LineString></Placemark>")
    s, e = blocks[si], blocks[ei]
    lines.append(
        f"<Placemark><name>START {cumd[si]:.3f} km</name>"
        f"<Point><coordinates>{s['lon']:.6f},{s['lat']:.6f},{s['alt']:.2f}</coordinates></Point></Placemark>"
    )
    lines.append(
        f"<Placemark><name>END {cumd[ei]:.3f} km</name>"
        f"<Point><coordinates>{e['lon']:.6f},{e['lat']:.6f},{e['alt']:.2f}</coordinates></Point></Placemark>"
    )
    lines.append("</Document></kml>")
    return "\n".join(lines)

@st.cache_data
def srt_zip(blocks, prefix):
    buf = io.BytesIO()
    z = zipfile.ZipFile(buf, 'w')
    cum = cumulative_dist([(b['lat'], b['lon']) for b in blocks])
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
    z.close()
    buf.seek(0)
    return buf.read()

# -- UI ------------------------------------------------------------------

# Step 1: Upload inputs
col1, col2 = st.columns(2)
with col1:
    up_srt = st.file_uploader("DJI .srt file", type='srt')
with col2:
    up_kml = st.file_uploader("Full-chainage KML", type=['kml', 'xml'])

if up_srt and up_kml:
    # Base chainage coordinates
    base_chain = parse_kml_2d(up_kml)
    if 'chain_rev' not in st.session_state:
        st.session_state.chain_rev = False

    # Apply reversal if needed
    chain_coords = list(reversed(base_chain)) if st.session_state.chain_rev else base_chain
    cum_chain = cumulative_dist(chain_coords)

    # Step 2: Generate & download markers
    if st.button("▶ Generate & Download 50 m Markers KML"):
        markers = generate_markers(chain_coords, cum_chain, 0.05)
        kml_m = kml_from_markers(markers)
        st.download_button("Download markers_50m.kml", kml_m, "markers_50m.kml",
                           mime="application/vnd.google-earth.kml+xml")

    # Step 3: Confirm direction
    dir_ok = st.radio("Are 50 m markers direction correct?", ["Yes", "Reverse Needed"])
    if dir_ok == "Reverse Needed":
        st.session_state.chain_rev = not st.session_state.chain_rev
        st.experimental_rerun()

    # Once confirmed:
    # Step 4: Parse SRT, project start/end
    blocks0 = parse_srt(up_srt)
    coords0 = [(b['lat'], b['lon']) for b in blocks0]
    alts0   = [b['alt'] for b in blocks0]
    cum0    = cumulative_dist(coords0)
    si0, start_km = project_to_line(coords0[0], chain_coords, cum_chain)
    ei0, end_km   = project_to_line(coords0[-1], chain_coords, cum_chain)

    # Shift SRT chainage so stationing starts at detected start
    offset = start_km
    cum0_shifted = [offset + d for d in cum0]

    st.write(f"Start → {start_km:.3f} km, End → {end_km:.3f} km, Total → {cum0[-1]*1000:.1f} m")

    # Step 5: Confirm SRT direction
    srt_dir = st.radio("Is SRT start→end correct?", ["Yes", "Reverse Needed"], key="srt_dir")
    if srt_dir == "Reverse Needed":
        blocks0.reverse()
        coords0.reverse()
        alts0.reverse()
        cum0 = cumulative_dist(coords0)
        cum0_shifted = [offset + d for d in cum0]
        si0, start_km = project_to_line(coords0[0], chain_coords, cum_chain)
        ei0, end_km   = project_to_line(coords0[-1], chain_coords, cum_chain)
        st.experimental_rerun()

    # Step 6: Visualize 3D
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=[c[1] for c in chain_coords],
        y=[c[0] for c in chain_coords],
        z=[0]*len(chain_coords),
        line=dict(color='gray'), name='Chainage Line'
    ))
    fig.add_trace(go.Scatter3d(
        x=[c[1] for c in coords0],
        y=[c[0] for c in coords0],
        z=alts0,
        line=dict(color='red'), name='SRT Route'
    ))
    fig.add_trace(go.Scatter3d(
        x=[coords0[0][1]],
        y=[coords0[0][0]],
        z=[alts0[0]],
        mode='markers+text',
        text=[f"Start {start_km:.3f} km"],
        marker=dict(size=5, color='green'),
        name='Start'
    ))
    fig.add_trace(go.Scatter3d(
        x=[coords0[-1][1]],
        y=[coords0[-1][0]],
        z=[alts0[-1]],
        mode='markers+text',
        text=[f"End {end_km:.3f} km"],
        marker=dict(size=5, color='blue'),
        name='End'
    ))
    fig.update_layout(
        scene=dict(xaxis_title='Lon', yaxis_title='Lat', zaxis_title='Alt (m)'),
        height=600
    )
    st.plotly_chart(fig)

    # Step 7: Export initial KML and SRT ZIP
    kml0 = kml_from_srt(blocks0, cum0_shifted, si0, ei0)
    st.download_button("Download initial_srt_route.kml", kml0, "initial_srt_route.kml",
                       mime="application/vnd.google-earth.kml+xml")
    zip0 = srt_zip(blocks0, "initial")
    st.download_button("Download initial SRT files ZIP", zip0, "initial_srt_files.zip",
                       mime="application/zip")

    # Step 8: Bulk processing
    up_bulk = st.file_uploader("Upload additional .srt files (bulk)", type='srt', accept_multiple_files=True)
    if up_bulk:
        all_buf = io.BytesIO()
        zall = zipfile.ZipFile(all_buf, 'w')
        for bf in up_bulk:
            name = bf.name.rsplit('.', 1)[0]
            blks = parse_srt(bf)
            crd = [(b['lat'], b['lon']) for b in blks]
            altb = [b['alt'] for b in blks]
            cumb = cumulative_dist(crd)
            si, sk = project_to_line(crd[0], chain_coords, cum_chain)
            ei, ek = project_to_line(crd[-1], chain_coords, cum_chain)
            offset_b = sk
            cumb_shift = [offset_b + d for d in cumb]
            kmlb = kml_from_srt(blks, cumb_shift, si, ei)
            zall.writestr(f"{name}_route.kml", kmlb)
            for fn in srt_zip(blks, name):  # reusing srt_zip to get filenames
                zall.write(fn)
        zall.close()
        all_buf.seek(0)
        st.download_button("Download all bulk outputs ZIP", all_buf.read(),
                           "all_bulk_outputs.zip", mime="application/zip")

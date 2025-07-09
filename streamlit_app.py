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
# Step 1: Upload & KML start‐chainage input
# ──────────────────────────────────────────────────────────────────────────

col1, col2, col3 = st.columns([2,2,1])
with col1:
    up_srt = st.file_uploader("DJI .srt file", type="srt")
with col2:
    up_kml = st.file_uploader("Full‐chainage KML", type=["kml","xml"])
with col3:
    chain_offset = st.number_input(
        "KML start chainage (km)", 
        min_value=0.0, step=0.001, format="%.3f", value=0.0
    )

if not (up_srt and up_kml):
    st.info("Please upload both an SRT file and your full‐chainage KML.")
    st.stop()

# ──────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────

def hav(p, q):
    R = 6371.0088
    φ1, φ2 = math.radians(p[0]), math.radians(q[0])
    dφ = math.radians(q[0] - p[0])
    dλ = math.radians(q[1] - p[1])
    a = math.sin(dφ/2)**2 + math.cos(φ1)*math.cos(φ2)*math.sin(dλ/2)**2
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))

@st.cache_data
def parse_srt(u):
    lines = u.getvalue().decode().splitlines()
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
def parse_kml_2d(u):
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    tree = ET.parse(io.BytesIO(u.getvalue()))
    coords = []
    for e in tree.findall('.//kml:LineString/kml:coordinates', ns):
        for part in e.text.strip().split():
            lon, lat, *_ = part.split(',')
            coords.append((float(lat), float(lon)))
    return coords

@st.cache_data
def cumulative_dist(coords, offset=0.0):
    cum = [offset]
    for A, B in zip(coords, coords[1:]):
        cum.append(cum[-1] + hav(A, B))
    return cum

@st.cache_data
def generate_markers(coords, cumd, interval_km=0.05):
    markers, d, idx = [], cumd[0], 0
    total = cumd[-1]
    while d <= total:
        while idx < len(cumd)-1 and cumd[idx+1] < d:
            idx += 1
        A, B = coords[idx], coords[idx+1]
        seg = cumd[idx+1] - cumd[idx]
        frac = (d - cumd[idx]) / seg if seg > 0 else 0
        markers.append((A[0] + frac*(B[0]-A[0]),
                        A[1] + frac*(B[1]-A[1]), d))
        d += interval_km
    return markers

def kml_header():
    return (
        "<?xml version='1.0' encoding='UTF-8'?>\n"
        "<!-- Author: Vijay Parmar; Community: BGol Community of Advanced Surveying and GIS Professionals -->\n"
        "<kml xmlns='http://www.opengis.net/kml/2.2'><Document>\n"
        "<Style id='chainStyle'>"
        "<LineStyle><color>ff0000ff</color><width>10.0</width></LineStyle>"
        "</Style>\n"
        "<Style id='srtStyle'>"
        "<LineStyle><color>ff0000ff</color><width>10.0</width></LineStyle>"
        "</Style>\n"
    )

@st.cache_data
def kml_from_markers(markers):
    xml = [kml_header()]
    for lat, lon, d in markers:
        xml.append(
            f"<Placemark><name>{d:.3f} km</name>"
            f"<Point><coordinates>{lon:.6f},{lat:.6f},0</coordinates></Point></Placemark>\n"
        )
    xml.append("</Document></kml>")
    return "".join(xml)

@st.cache_data
def project_to_line(pt, coords, cumd):
    dists = [hav(pt, c) for c in coords]
    i = min(range(len(dists)), key=lambda i: dists[i])
    return i, cumd[i]

# …—— PATCHED kml_from_srt ————————————————
@st.cache_data
def kml_from_srt(blocks, start_km, end_km):
    xml = [kml_header()]
    # full SRT LineString
    xml.append(
        "<Placemark><styleUrl>#srtStyle</styleUrl><name>SRT Route</name>"
        "<LineString><coordinates>\n"
    )
    for b in blocks:
        xml.append(f"{b['lon']:.6f},{b['lat']:.6f},{b['alt']:.2f}\n")
    xml.append("</coordinates></LineString></Placemark>\n")
    # START marker at first block
    first = blocks[0]
    xml.append(
        f"<Placemark><styleUrl>#srtStyle</styleUrl>"
        f"<name>START {start_km:.3f} km</name>"
        f"<Point><coordinates>"
        f"{first['lon']:.6f},{first['lat']:.6f},{first['alt']:.2f}"
        f"</coordinates></Point></Placemark>\n"
    )
    # END marker at last block
    last = blocks[-1]
    xml.append(
        f"<Placemark><styleUrl>#srtStyle</styleUrl>"
        f"<name>END {end_km:.3f} km</name>"
        f"<Point><coordinates>"
        f"{last['lon']:.6f},{last['lat']:.6f},{last['alt']:.2f}"
        f"</coordinates></Point></Placemark>\n"
    )
    xml.append("</Document></kml>")
    return "".join(xml)
# …———————————————————————————————————————

def srt_mapping(blocks, prefix, offset=0.0):
    cum = cumulative_dist([(b['lat'],b['lon']) for b in blocks], offset)
    labels = [
        "01_Latitude","02_Longitude","03_Altitude",
        "04_Timer","05_LatLonOutput","06_Chainage"
    ]
    cols = [
        [f"{b['lat']:.6f}" for b in blocks],
        [f"{b['lon']:.6f}" for b in blocks],
        [f"{b['alt']:.2f}" for b in blocks],
        [b['tim'] for b in blocks],
        [f"{b['lat']:.6f}, {b['lon']:.6f}" for b in blocks],
        [f"{d:.3f} km" for d in cum]
    ]
    mapping = {}
    for lab, col in zip(labels, cols):
        fn = f"{prefix}_{lab}.srt"
        txt = "\n\n".join(
            f"{blocks[i]['idx']}\n{blocks[i]['range']}\n{col[i]}"
            for i in range(len(blocks))
        )
        mapping[fn] = txt
    return mapping

# ──────────────────────────────────────────────────────────────────────────
# Step 2: Markers + FLIP + Proceed
# ──────────────────────────────────────────────────────────────────────────

st.subheader("Step 2: Generate 50 m Markers")
flip = st.checkbox("FLIP chainage direction", key="flip_chain")
base_chain = parse_kml_2d(up_kml)
chain_coords = list(reversed(base_chain)) if flip else base_chain
cum_chain    = cumulative_dist(chain_coords, offset=chain_offset)

if st.button("▶ Generate & Download 50 m Markers KML"):
    kmlm = kml_from_markers(generate_markers(chain_coords, cum_chain))
    st.download_button("Download markers_50m.kml", kmlm,
                       "markers_50m.kml","application/vnd.google-earth.kml+xml")

if "proceed" not in st.session_state:
    st.session_state.proceed = False
if st.button("Proceed further"):
    st.session_state.proceed = True

# ──────────────────────────────────────────────────────────────────────────
# Steps 4–9 after proceed
# ──────────────────────────────────────────────────────────────────────────

if st.session_state.proceed:
    # Step 4: Parse & project SRT
    blocks0 = parse_srt(up_srt)
    pts0    = [(b['lat'],b['lon']) for b in blocks0]
    alts0   = [b['alt'] for b in blocks0]
    cum0    = cumulative_dist(pts0)
    si0, sk0 = project_to_line(pts0[0],   chain_coords, cum_chain)
    ei0, ek0 = project_to_line(pts0[-1],  chain_coords, cum_chain)

    st.subheader("Step 4: SRT Stationing")
    st.write(f"• Start → **{sk0:.3f} km**")
    st.write(f"• End   → **{ek0:.3f} km**")
    st.write(f"• Total → **{cum0[-1]*1000:.1f} m**")

    # Step 6: 3D Visualization
    st.subheader("Step 6: 3D Visualization")
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=[c[1] for c in chain_coords],
        y=[c[0] for c in chain_coords],
        z=[0]*len(chain_coords),
        line=dict(color='blue', width=10.0), name='Chainage'
    ))
    fig.add_trace(go.Scatter3d(
        x=[p[1] for p in pts0],
        y=[p[0] for p in pts0],
        z=alts0,
        line=dict(color='red', width=4.0), name='SRT Route'
    ))
    fig.add_trace(go.Scatter3d(
        x=[pts0[0][1]], y=[pts0[0][0]], z=[alts0[0]],
        mode='markers+text',
        text=[f"Start {sk0:.3f} km"],
        marker=dict(size=5, color='green'), name='Start'
    ))
    fig.add_trace(go.Scatter3d(
        x=[pts0[-1][1]], y=[pts0[-1][0]], z=[alts0[-1]],
        mode='markers+text',
        text=[f"End {ek0:.3f} km"],
        marker=dict(size=5, color='blue'), name='End'
    ))
    fig.update_layout(
        scene=dict(xaxis_title='Lon', yaxis_title='Lat', zaxis_title='Alt (m)'),
        height=600
    )
    st.plotly_chart(fig)

    # Step 7: Export initial KML & SRT ZIP
    st.subheader("Step 7: Export Initial Files")
    kml0 = kml_from_srt(blocks0, sk0, ek0)
    st.download_button("Download initial_srt_route.kml", kml0,
                       "initial_srt_route.kml","application/vnd.google-earth.kml+xml")

    map0 = srt_mapping(blocks0, "initial", offset=sk0)
    buf0 = io.BytesIO(); z0 = zipfile.ZipFile(buf0, 'w')
    for fn, txt in map0.items():
        z0.writestr(fn, txt)
    z0.close(); buf0.seek(0)
    st.download_button("Download initial SRT files ZIP", buf0.read(),
                       "initial_srt_files.zip","application/zip")

    # Step 8: Bulk Processing
    st.subheader("Step 8: Bulk Processing")
    up_bulk = st.file_uploader(
        "Upload additional .srt files (bulk)",
        type='srt', accept_multiple_files=True
    )
    if up_bulk:
        all_buf = io.BytesIO(); zall = zipfile.ZipFile(all_buf,'w')
        for bf in up_bulk:
            name = bf.name.rsplit('.',1)[0]
            blks2 = parse_srt(bf)
            pts2  = [(b['lat'],b['lon']) for b in blks2]
            cum2  = cumulative_dist(pts2)
            s2, sk2 = project_to_line(pts2[0],  chain_coords, cum_chain)
            e2, ek2 = project_to_line(pts2[-1], chain_coords, cum_chain)
            kml2 = kml_from_srt(blks2, sk2, ek2)
            zall.writestr(f"{name}_route.kml", kml2)
            map2 = srt_mapping(blks2, name, offset=sk2)
            for fn, txt in map2.items():
                zall.writestr(fn, txt)
        zall.close(); all_buf.seek(0)
        st.download_button("Download all bulk outputs ZIP", all_buf.read(),
                           "all_bulk_outputs.zip","application/zip")

    st.success("✅ All outputs are ready. - BGol")

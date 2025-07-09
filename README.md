# Drone Chainage & SRT Processing Streamlit App

A Streamlit application that automates extraction and processing of GPS metadata from DJI drone video `.srt` files, generates chainage markers along a reference alignment (KML), and produces synchronized SRT outputs and 3D visualizations.

---

## 🚀 Purpose

* **Extract** latitude, longitude, altitude, and timestamp data embedded in DJI `.srt` files.
* **Chainage**: compute distances along a surveyed alignment (a 2D KML LineString) and insert regular markers (every 50 m).
* **Project** drone track start/end onto the chainage alignment to find matching station values.
* **Visualize** the drone path and alignment in an interactive 3D view.
* **Export**:

  * 3D KMLs of the drone route with START/END placemarks.
  * Six separate SRT files (Latitude, Longitude, Altitude, Timer, LatLon, Chainage).
  * Bulk‐processing support for multiple `.srt` files in one go.

---

## 📋 Features & Workflow

1. **Upload Inputs**

   * DJI `.srt` (with `[latitude: …]`, `[longitude: …]`, `[altitude: …]`).
   * Full‐chainage `.kml` (2D LineString of reference alignment).

2. **Generate 50 m Markers**

   * Interpolate along the KML every 0.05 km.
   * Download `markers_50m.kml` for inspection.

3. **Confirm Markers Direction**

   * If backward, reverse the alignment and regenerate markers.

4. **Project SRT Start/End**

   * Parse the initial `.srt` into GPS track points.
   * Compute total route length.
   * Find nearest chainage-station for the first/last points.

5. **Confirm SRT Direction**

   * Reverse the track if start/end mapping is inverted.

6. **Visualize in 3D**

   * Chainage line (gray) at altitude = 0.
   * Drone track (red) in true altitude.
   * START (green) and END (blue) markers with station labels.

7. **Export Initial SRT-Derived Files**

   * Download `initial_srt_route.kml` with START/END.
   * Download a ZIP of six SRT files prefixed `initial_01_…`–`initial_06_…`.

8. **Bulk Processing** (optional)

   * Upload multiple additional `.srt` files.
   * Skip Steps 2–6 and generate 3D KML + SRT ZIP for each.

9. **Bundle & Download**

   * All KMLs and SRTs are packaged into a final `all_outputs.zip` for one-click retrieval.

---

## ⚙️ Installation & Running Locally

1. **Clone the repo**:

   ```bash
   git clone https://github.com/<your-user>/drone-chainage-app.git
   cd drone-chainage-app
   ```

2. **Install dependencies**:

   ```bash
   pip install -r requirements.txt
   ```

   *(Requires: Python 3.7+, streamlit, plotly)*

3. **Run**:

   ```bash
   streamlit run streamlit_app.py
   ```

4. **Access** the app in your browser at `http://localhost:8501`.

---

## ☁️ Deploying on Streamlit Cloud

1. Push this repo to GitHub.
2. Go to [streamlit.io/cloud](https://streamlit.io/cloud), connect your GitHub account.
3. Create a new app, select this repo and `streamlit_app.py` as the entrypoint.
4. Click **Deploy**.

---

## 🔍 Under the Hood: Mathematical Background

### 1. Haversine Formula

Computes great‐circle distance between two GPS points $(φ_1, λ_1)$ and $(φ_2, λ_2)$:

$$
Δφ = φ_2 - φ_1,
\quad Δλ = λ_2 - λ_1,
\quad a = \sin^2(Δφ/2) + \cos φ_1 \, \cos φ_2 \, \sin^2(Δλ/2),
\quad d = 2R \arcsin(\sqrt{a}),
$$

where:

* $φ$ = latitude (radians), $λ$ = longitude (radians)
* $R$ ≈ 6371 km (Earth’s radius)
* $d$ = distance along Earth’s surface.

### 2. Cumulative Chainage

Given points $P_0,P_1,…,P_n$:

$$
C_k = \sum_{i=1}^k d(P_{i-1},P_i),
$$

where $d()$ is the segment length from Haversine.

### 3. Regular Interval Markers

To place markers every $Δ$ km:

1. Build cumulative array $[C_0=0,C_1,…,C_n]$.
2. For each$ D=0, Δ, 2Δ,…,C_n$: find segment $i$ with
   $C_{i-1} < D ≤ C_i$, compute fraction
   $f = (D - C_{i-1})/(C_i - C_{i-1})$, and
   interpolated point $P(D) = (1-f)P_{i-1} + fP_i$.

---

## 📜 License & Credits

* **Author**: Vijay Parmar (@VIJAYPARMAR)
* **License**: Open Source 
* **Community**: BGol Community of Advanced Surveying and GIS Professionals

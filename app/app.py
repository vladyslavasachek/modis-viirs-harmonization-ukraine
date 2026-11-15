"""Harmonized burning-activity calendar for Ukraine (MODIS + VIIRS) — Streamlit app.
Run: streamlit run app/app.py   (uses only data/processed/)"""
import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import shapely
import streamlit as st

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import lib  # noqa: E402
from lib import PNG_CONFIG, UNIT  # noqa: E402

st.set_page_config(page_title="Ukraine burning calendar", page_icon="🔥", layout="wide")


# ============================================================================ sidebar: AOI, domain, dates
def sidebar():
    st.sidebar.title("Area & period")
    mode = st.sidebar.radio("Area of interest (AOI)", ["Preset", "Point + radius", "Draw polygon"], horizontal=False)
    aoi = {}
    if mode == "Preset":
        presets = lib.PRESETS_FIXED[:3] + lib.oblast_names()
        name = st.sidebar.selectbox("Preset", presets, index=0)
        aoi = dict(kind="preset", name=name, cells=lib.aoi_cells(name))
    elif mode == "Point + radius":
        c1, c2 = st.sidebar.columns(2)
        lat = c1.number_input("Latitude", 44.0, 53.0, 50.45, 0.05, format="%.3f")
        lon = c2.number_input("Longitude", 22.0, 41.0, 30.52, 0.05, format="%.3f")
        r = st.sidebar.slider("Radius (km)", 5, 400, 100, 5)
        aoi = dict(kind="circle", lat=lat, lon=lon, r=r, cells=lib.circle_cells(lat, lon, r),
                   name=f"{r} km around {lat:.2f}N {lon:.2f}E")
    else:
        poly = st.session_state.get("drawn_polygon")
        aoi = dict(kind="polygon", polygon=poly, cells=lib.aoi_cells("Custom", poly) if poly is not None else np.array([]),
                   name="Custom polygon")
        st.sidebar.info("Draw a polygon on the map on the Calendar page.")
    domain_label = st.sidebar.radio("Spatial domain", list(lib.DOMAIN_LABELS.values()), index=0,
                                    help=lib.CONFLICT_NOTE)
    domain = [k for k, v in lib.DOMAIN_LABELS.items() if v == domain_label][0]
    weeks = lib.all_weeks()
    dmin, dmax = weeks["week_start"].min().date(), lib.load_detections()["datetime_utc"].max().date()
    dr = st.sidebar.date_input("Date range", (pd.Timestamp("2001-01-01").date(), dmax), min_value=dmin, max_value=dmax)
    if not isinstance(dr, (tuple, list)) or len(dr) != 2:
        dr = (dmin, dmax)
    small = len(aoi["cells"]) < 20
    base_opts = ["2001–2025 (incl. MODIS-based estimates before 2012)", "2012–2025 (observed VIIRS only)"]
    bsel = st.sidebar.radio("Baseline years (anomalies, KPI)", base_opts, index=1 if small else 0,
                            key=f"baseline_{'small' if small else 'large'}",
                            help="For small AOIs (<20 cells of 0.1°) the default is VIIRS-only: before 2012 the "
                                 "MODIS-based estimate is zero wherever MODIS saw nothing, which distorts the "
                                 "distribution of weekly values in small areas.")
    base_years = list(range(2012, 2026)) if bsel.startswith("2012") else list(range(2001, 2026))
    st.sidebar.caption(f"AOI: **{aoi['name']}** · {len(aoi['cells'])} cells of 0.1°")
    st.sidebar.caption(f"Domain: {domain_label}")
    return aoi, domain, pd.Timestamp(dr[0]), pd.Timestamp(dr[1]), base_years


def get_series(aoi, domain, d0, d1):
    if len(aoi["cells"]) == 0:
        return None
    wk = lib.aoi_weekly(tuple(sorted(map(int, aoi["cells"]))), domain)
    scale = wk.attrs.get("scale")
    wk = wk[(wk["week_start"] >= d0 - pd.Timedelta(days=6)) & (wk["week_start"] <= d1)].copy()
    wk.attrs["scale"] = scale
    return wk


def hover_text(df):
    return lib.fmt_interval(df) + "<br>" + df.get("source_label", pd.Series("", index=df.index)).astype(str)


# ============================================================================ pages
def page_calendar():
    aoi, domain, d0, d1 = AOI, DOMAIN, D0, D1
    st.title("Burning-activity calendar")
    st.caption(f"{aoi['name']} · {lib.DOMAIN_LABELS[domain]} · unit: **{UNIT}** (S-NPP 375 m scale). "
               "Detections are satellite observations, not a count of fires.")
    if aoi["kind"] == "polygon":
        draw_map()
    if len(aoi["cells"]) == 0:
        st.warning("The AOI contains no 0.1° cells with detections (or no polygon drawn yet).")
        return
    wk = get_series(aoi, domain, d0, d1)
    det = lib.aoi_detections(aoi, domain)
    calib_ta, calib_t, n20q, _ = lib.load_calib()

    # ---- KPI: current month-to-date vs historical median
    k = lib.month_to_date_kpi(det, n20q, base_years=BASE)
    by = f"{BASE[0]}–{BASE[-1]}"
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(f"This month to date ({k['start']:%d %b}–{k['end']:%d %b %Y})", f"{k['value']:,.0f}",
              help=f"95% interval {k['value_lo']:,.0f}–{k['value_hi']:,.0f} (NOAA-20→S-NPP calibration). NRT data.")
    c2.metric(f"Historical median, same window ({by})", f"{k['median']:,.0f}",
              help=f"95% bootstrap interval of the median: {k['median_lo']:,.0f}–{k['median_hi']:,.0f}. "
                   "Years before 2012 are MODIS-based estimates.")
    c3.metric(f"Percentile vs {by}", f"{k['percentile']:.0f}th")
    c4.metric("Ratio to median", f"{k['value'] / k['median']:.2f}×" if k["median"] > 0 else "n/a")
    st.caption(f"Same calendar window ({k['start']:%d %b}–{k['end']:%d %b}) in every baseline year {by}. "
               f"Current value 95% interval {k['value_lo']:,.0f}–{k['value_hi']:,.0f}; median 95% interval "
               f"{k['median_lo']:,.0f}–{k['median_hi']:,.0f}. **The current month is near-real-time (NRT) data** "
               "(NOAA-20 calibrated to S-NPP after S-NPP ended on 2 Nov 2026); NRT can differ slightly from the "
               "later standard archive. Baseline years before 2012 are MODIS-based estimates.")

    # ---- weekly series with intervals
    fig = go.Figure()
    x = wk["week_start"]
    fig.add_trace(go.Scatter(x=x, y=wk["hi95"], mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=x, y=wk["lo95"], fill="tonexty", mode="lines", line=dict(width=0), name="95% interval",
                             fillcolor="rgba(33,145,140,0.2)", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=x, y=wk["hi80"], mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=x, y=wk["lo80"], fill="tonexty", mode="lines", line=dict(width=0), name="80% interval",
                             fillcolor="rgba(33,145,140,0.4)", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=x, y=wk["point"], line=dict(color="#21918c", width=1.4), name="Harmonized",
                             text=hover_text(wk), hovertemplate="%{x|%d %b %Y}: %{text}<extra></extra>"))
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10), yaxis_title=f"{UNIT} / week",
                      title="Weekly harmonized activity", legend=dict(orientation="h", y=1.12))
    st.plotly_chart(fig, use_container_width=True, config=PNG_CONFIG)

    # ---- heatmaps
    c1 = c2 = st.container()
    hm = wk.pivot_table(index="year", columns="woy", values="point", aggfunc="sum")
    txt = wk.assign(t=hover_text(wk)).pivot_table(index="year", columns="woy", values="t", aggfunc="first")
    f1 = go.Figure(go.Heatmap(z=hm.values, x=hm.columns + 1, y=hm.index, colorscale="Viridis",
                              text=txt.reindex_like(hm).values,
                              hovertemplate="year %{y}, week %{x}<br>%{text}<extra></extra>",
                              colorbar=dict(title="hot spots")))
    f1.update_layout(height=480, title="Week of year × year", xaxis_title="week of year", yaxis_title="year",
                     margin=dict(l=10, r=10, t=40, b=10), yaxis=dict(autorange="reversed"))
    c1.plotly_chart(f1, use_container_width=True, config=PNG_CONFIG)
    mon = lib.aggregate(wk, ["year", "month"])
    hm2 = mon.pivot(index="year", columns="month", values="point")
    txt2 = mon.assign(t=lib.fmt_interval(mon)).pivot(index="year", columns="month", values="t")
    f2 = go.Figure(go.Heatmap(z=hm2.values, x=[pd.Timestamp(2000, m, 1).strftime("%b") for m in hm2.columns],
                              y=hm2.index, colorscale="Viridis", text=txt2.values,
                              hovertemplate="%{x} %{y}<br>%{text}<extra></extra>",
                              colorbar=dict(title="hot spots")))
    f2.update_layout(height=480, title="Month × year", yaxis=dict(autorange="reversed"),
                     margin=dict(l=10, r=10, t=40, b=10))
    c2.plotly_chart(f2, use_container_width=True, config=PNG_CONFIG)
    st.caption("Hover shows the value with 80/95% intervals (MODIS-based years) or 'observed' (S-NPP years) and "
               "the data source. Use the camera icon on any plot to download it as PNG.")

    # ---- annual series + trend
    yr = lib.aggregate(wk, ["year"])
    full = yr[(yr["year"] >= max(2001, d0.year + (d0.dayofyear > 7))) & (yr["year"] <= min(2025, d1.year - (d1.dayofyear < 358)))]
    fa = go.Figure()
    fa.add_trace(go.Bar(x=yr["year"], y=yr["point"], marker_color="#31688e", name="Harmonized",
                        error_y=dict(type="data", symmetric=False, array=yr["hi95"] - yr["point"],
                                     arrayminus=yr["point"] - yr["lo95"], color="#555"),
                        hovertext=lib.fmt_interval(yr), hovertemplate="%{x}: %{hovertext}<extra></extra>"))
    det_all = lib.aoi_detections(aoi, "all")
    # trend on a constant area: "excl" -> conflict-affected oblasts removed in all years (north 2022 window kept)
    tr_ann = lib.annual_from_weights(lib.constant_area(det_all) if domain == "excl" else det_all)["harmonized"]
    tr_years = full["year"].astype(int).values
    if len(full) >= 8:
        t = lib.sen_mk(tr_years, tr_ann.reindex(tr_years, fill_value=0).values)
        fa.add_trace(go.Scatter(x=tr_years, y=t["intercept"] + t["slope"] * tr_years, mode="lines",
                                line=dict(color="#fde725", width=3, dash="dash"), name="Sen's slope"))
        st.markdown(f"**Long-term trend ({tr_years.min()}–{tr_years.max()}, complete years, constant area):** "
                    f"Sen's slope {t['slope']:+,.0f} {UNIT}/year (95% CI {t['lo']:+,.0f} to {t['hi']:+,.0f}; "
                    f"{t['rel']:+.1f}%/yr of the median); Mann–Kendall τ = {t['tau']:.2f}, p = {t['p']:.3f}. "
                    "Descriptive only: no correction for autocorrelation, and pre-2012 values are model estimates "
                    "(their intervals are not propagated into the slope CI).")
    with st.expander("Trend robustness: harmonized vs single-sensor raw series", expanded=True):
        rob = lib.trend_robustness(det_all)
        signs = np.sign(rob["Sen's slope / yr"].iloc[:3].values)
        st.dataframe(rob.style.format({"Sen's slope / yr": "{:+,.0f}", "% of median / yr": "{:+.1f}%",
                                       "Mann–Kendall τ": "{:.2f}", "p": "{:.3f}"}),
                     hide_index=True, use_container_width=True)
        if len(set(signs[signs != 0])) > 1:
            st.warning("The three series disagree in the sign of the trend; treat the harmonized trend as not robust.")
        else:
            st.caption("All three constant-area series agree in sign. Raw single-sensor slopes are in that sensor's "
                       "own detection units; compare signs and % of median per year, not absolute values. Fixed "
                       "periods, complete years, independent of the date-range filter. Compare the 2001–2025 and "
                       "pre-invasion (2001–2021) rows: a trend that is significant only with 2022–2025 included "
                       "is driven by the war-time years.")
    fa.update_layout(height=330, title="Annual harmonized activity (error bars: 95% interval)",
                     yaxis_title=f"{UNIT} / year", margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fa, use_container_width=True, config=PNG_CONFIG)

    # ---- top months / years
    c1, c2 = st.columns(2)
    tm = mon.sort_values("point", ascending=False).head(10)
    c1.markdown("**Top 10 months**")
    c1.dataframe(pd.DataFrame({"month": [f"{pd.Timestamp(2000, int(m), 1):%b} {int(y)}" for y, m in zip(tm["year"], tm["month"])],
                               UNIT: lib.fmt_interval(tm).values}), hide_index=True, use_container_width=True)
    ty = yr.sort_values("point", ascending=False).head(10)
    c2.markdown("**Top 10 years**")
    c2.dataframe(pd.DataFrame({"year": ty["year"].astype(int).values, UNIT: lib.fmt_interval(ty).values}),
                 hide_index=True, use_container_width=True)
    st.caption("Partial years (2000, 2026) and partial date ranges are included as-is.")

    # ---- downloads
    c1, c2 = st.columns(2)
    cal = wk[["week_start", "year", "woy", "month", "point", "lo80", "hi80", "lo95", "hi95", "source_label",
              "raw_modis", "raw_snpp", "raw_noaa20"]].rename(columns={"point": "viirs_equivalent_hotspots",
                                                                      "woy": "week_of_year", "source_label": "source"})
    c1.download_button("⬇ Calendar CSV (weekly, with intervals)", cal.to_csv(index=False).encode(),
                       "calendar_weekly.csv", "text/csv")
    dsel = det[(det["datetime_utc"] >= d0) & (det["datetime_utc"] < d1 + pd.Timedelta(days=1))]
    key = (aoi["name"], domain, str(d0), str(d1), len(dsel))
    if st.session_state.get("hs_key") != key:
        st.session_state.pop("hs_csv", None)
    if "hs_csv" not in st.session_state:
        if c2.button(f"Prepare hot spots CSV ({len(dsel):,} detections in AOI & date range)"):
            st.session_state["hs_csv"] = dsel[["sensor", "datetime_utc", "lat", "lon", "frp", "confidence", "daynight",
                                               "processing", "conflict_flag", "oblast"]].to_csv(index=False).encode()
            st.session_state["hs_key"] = key
            st.rerun()
    else:
        c2.download_button(f"⬇ Hot spots CSV ({len(dsel):,} detections)", st.session_state["hs_csv"],
                           "hotspots.csv", "text/csv")


def draw_map():
    import folium
    from folium.plugins import Draw
    from streamlit_folium import st_folium
    m = folium.Map(location=[49.0, 31.0], zoom_start=6, tiles="OpenStreetMap")
    folium.GeoJson(lib.load_boundaries(), style_function=lambda f: dict(
        color="#21918c" if f["properties"].get("kind") == "oblast" else "#e66101", weight=1, fillOpacity=0),
        tooltip=folium.GeoJsonTooltip(fields=["shapeName"], labels=False)).add_to(m)
    Draw(export=False, draw_options=dict(polyline=False, circle=False, marker=False, circlemarker=False,
                                         rectangle=True, polygon=True)).add_to(m)
    out = st_folium(m, height=380, use_container_width=True, returned_objects=["last_active_drawing"])
    g = (out or {}).get("last_active_drawing")
    if g and g.get("geometry"):
        poly = shapely.geometry.shape(g["geometry"])
        if st.session_state.get("drawn_wkt") != poly.wkt:
            st.session_state["drawn_polygon"] = poly
            st.session_state["drawn_wkt"] = poly.wkt
            st.rerun()


def page_anomalies():
    import calendar_stats as cs
    aoi, domain = AOI, DOMAIN
    st.title("Anomalies & critical periods")
    st.caption(f"{aoi['name']} · {lib.DOMAIN_LABELS[domain]} · unit: **{UNIT}**. Baseline: complete years "
               f"{BASE[0]}–{BASE[-1]} (sidebar), always excluding the year being evaluated. Anomalies are percentiles, not z-scores "
               "(the data are skewed and zero-inflated). Descriptive statistics only — no forecasts.")
    if len(aoi["cells"]) == 0:
        st.warning("Empty AOI.")
        return
    wk = lib.aoi_weekly(tuple(sorted(map(int, aoi["cells"]))), domain)
    mon = lib.aggregate(wk, ["year", "month"])
    years = sorted(wk["year"].unique())
    ysel = st.selectbox("Year to evaluate", years[::-1], index=years[::-1].index(2025) if 2025 in years else 0)

    # ---- climatology with selected year
    cl = cs.climatology(wk, "woy", exclude_year=ysel, years=BASE)
    cur = wk[wk["year"] == ysel].set_index("woy").reindex(range(52))
    x = np.arange(1, 53)
    f = go.Figure()
    f.add_trace(go.Scatter(x=x, y=cl["p90"], mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip"))
    f.add_trace(go.Scatter(x=x, y=cl["p10"], mode="lines", line=dict(width=0), fill="tonexty",
                           fillcolor="rgba(120,120,120,0.3)", name="10–90% range (baseline)"))
    f.add_trace(go.Scatter(x=x, y=cl["median"], mode="lines", line=dict(color="#808080", dash="dash"),
                           name="Baseline median"))
    est = cur["share_estimated"].fillna(0) > 0.001
    f.add_trace(go.Scatter(x=x, y=cur["point"], mode="lines+markers", line=dict(color="#21918c", width=2.5),
                           marker=dict(size=5), name=str(ysel),
                           error_y=dict(type="data", symmetric=False, visible=bool(est.any()),
                                        array=(cur["hi95"] - cur["point"]).where(est, 0),
                                        arrayminus=(cur["point"] - cur["lo95"]).where(est, 0)),
                           text=lib.fmt_interval(cur.fillna(0)), hovertemplate="week %{x}: %{text}<extra></extra>"))
    f.update_layout(height=360, title=f"Weekly activity in {ysel} vs climatology of the other years",
                    xaxis_title="week of year", yaxis_title=f"{UNIT} / week", margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(f, use_container_width=True, config=PNG_CONFIG)

    # ---- percentile anomalies
    st.subheader("Percentile anomalies")
    gran = st.radio("Resolution", ["Month", "Week"], horizontal=True)
    period = "month" if gran == "Month" else "woy"
    src = mon if period == "month" else wk
    P, V = cs.percentile_anomalies(src, period, years=BASE)
    xs = ([pd.Timestamp(2000, m, 1).strftime("%b") for m in P.columns] if period == "month"
          else [c + 1 for c in P.columns])
    txt = np.vectorize(lambda v, p: f"{v:,.0f} hot spots · {p:.0f}th percentile")(V.values, P.values)
    fh = go.Figure(go.Heatmap(z=P.values, x=xs, y=P.index, colorscale="RdBu_r", zmid=50, zmin=0, zmax=100,
                              text=txt, hovertemplate="%{y} %{x}: %{text}<extra></extra>",
                              colorbar=dict(title="percentile")))
    fh.update_layout(height=520, title=f"{gran}ly activity as percentile of the same {gran.lower()} in other years "
                     "(blue = below, red = above median)", yaxis=dict(autorange="reversed"),
                     margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fh, use_container_width=True, config=PNG_CONFIG)
    st.caption(f"Percentiles are computed against complete baseline years {BASE[0]}–{BASE[-1]} without the evaluated year; for "
               "2000 and 2026 (partial) only the covered weeks are meaningful. Pre-2012 values are MODIS-based "
               "estimates; their intervals are shown in the Calendar page hover.")

    # ---- critical periods
    st.subheader("Critical periods")
    cp, thr = cs.critical_periods(wk, years=BASE)
    fc = go.Figure(go.Bar(x=cp["woy"] + 1, y=cp["prob"] * 100,
                          marker_color=np.where(cp["critical"], "#b2182b", "#92c5de"),
                          error_y=dict(type="data", symmetric=False, array=(cp["hi95"] - cp["prob"]) * 100,
                                       arrayminus=(cp["prob"] - cp["lo95"]) * 100, color="#555"),
                          hovertemplate="week %{x}: %{y:.0f}%<extra></extra>"))
    fc.add_hline(y=50, line_dash="dash", line_color="#444", annotation_text="50%")
    fc.update_layout(height=330, xaxis_title="week of year", yaxis_title="probability (%)",
                     title=f"Probability that a week exceeds the baseline 75th percentile of weekly activity "
                           f"(≈{thr:,.0f} {UNIT}/week); red = critical (>50%)", margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fc, use_container_width=True, config=PNG_CONFIG)
    crit = cp.loc[cp["critical"], "woy"].values
    if len(crit):
        dates = [(pd.Timestamp(2001, 1, 1) + pd.Timedelta(days=int(w) * 7)).strftime("%d %b") for w in crit]
        st.markdown("**Critical weeks:** " + ", ".join(f"week {w + 1} (from {d})" for w, d in zip(crit, dates)))
    st.caption("The 75th-percentile threshold is pooled over all weeks of the baseline years (computed without the "
               "evaluated year); a per-week threshold would be exceeded in ~25% of years by construction. "
               "95% intervals: bootstrap over years. For sparse AOIs the threshold can be 0, i.e. 'any activity'.")

    # ---- season metrics
    st.subheader("Season timing")
    split = cs.season_split(wk, years=BASE)
    sm = cs.season_metrics(wk[wk["year"].between(2001, 2025)], split=split)
    st.caption(("Bimodal season detected: spring / summer-autumn split at week "
                f"{split + 1} (~{(pd.Timestamp(2001, 1, 1) + pd.Timedelta(days=split * 7)):%d %b}).")
               if split is not None else "No clear bimodal season; metrics computed over the whole year.")
    fs = go.Figure()
    colors = {"spring": "#1b9e77", "summer-autumn": "#d95f02", "year": "#7570b3"}
    for sea, g in sm.groupby("season"):
        c = colors.get(sea, "#333")
        fs.add_trace(go.Scatter(x=g["year"], y=g["peak_doy"], mode="markers", marker=dict(color=c, size=8),
                                name=f"{sea}: peak",
                                error_y=dict(type="data", symmetric=False, array=g["end_doy"] - g["peak_doy"],
                                             arrayminus=g["peak_doy"] - g["onset_doy"], color=c),
                                hovertemplate="%{x}: peak day %{y}<extra></extra>"))
    fs.update_layout(height=380, yaxis_title="day of year (bars: onset 10% → end 90%)",
                     title="Season onset, peak and end by year", margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fs, use_container_width=True, config=PNG_CONFIG)
    st.dataframe(sm.groupby("season")[["onset_doy", "peak_doy", "end_doy", "length_days"]].median()
                 .rename(columns=lambda c: "median " + c.replace("_", " ")), use_container_width=True)
    c1, c2 = st.columns(2)
    c1.download_button("⬇ Percentile anomalies CSV", P.round(1).to_csv().encode(), f"anomaly_percentiles_{period}.csv")
    c2.download_button("⬇ Season metrics CSV", sm.to_csv(index=False).encode(), "season_metrics.csv")


def page_harmonization():
    st.title("Harmonization: MODIS → VIIRS S-NPP equivalent")
    ann = pd.read_csv(os.path.join(lib.PROC, "ukraine_annual_harmonized.csv"))
    a = ann[(ann["year"] >= 2001) & (ann["year"] <= 2025)]
    fig = go.Figure()
    m = a[a["year"] <= 2011]
    fig.add_trace(go.Scatter(x=m["year"], y=m["hi95"], mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=m["year"], y=m["lo95"], fill="tonexty", mode="lines", line=dict(width=0), name="95% interval",
                             fillcolor="rgba(49,104,142,0.2)"))
    fig.add_trace(go.Scatter(x=a["year"], y=a["Terra"] + a["Aqua"], name="Raw MODIS (Terra+Aqua)",
                             line=dict(color="#666666")))
    s = a[a["year"] >= 2012]
    fig.add_trace(go.Scatter(x=s["year"], y=s["SNPP"], name="Raw VIIRS S-NPP", line=dict(color="#e66101")))
    fig.add_trace(go.Scatter(x=a["year"], y=a["point"], name="Harmonized", line=dict(color="#31688e", width=3)))
    fig.update_layout(height=420, yaxis_title=f"{UNIT} / year (raw lines: detections)",
                      title="All of Ukraine: the 2012 sensor jump disappears after harmonization",
                      margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fig, use_container_width=True, config=PNG_CONFIG)
    st.markdown("""
**Method.** Reference scale = VIIRS S-NPP 375 m. For 2000–2012 (and S-NPP outage weeks before NOAA-20) the
S-NPP-equivalent value is estimated from MODIS detections with **B1**: a ratio VIIRS/MODIS by month × day/night,
fitted on the overlap 2012–2021 (pre-invasion), shrunk towards the global ratio for small strata. Nov 2000 –
early Jul 2002 (Terra only) uses a separate Terra-only ratio and is flagged *lower confidence*. After S-NPP ends
(2 Nov 2026) and in S-NPP outage weeks from 2018, NOAA-20 is used × 1.011 (calibrated to S-NPP).

**Model choice by leave-one-year-out validation** (2012–2021, predicting S-NPP from MODIS only). B0 is slightly
better at the yearly level, but B1 wins at the month level, which is what the calendar uses. The Negative-Binomial
GLM (M1) did not validate better than the simple ratios.""")
    v = pd.concat([pd.read_csv(os.path.join(lib.PROC, "validation_loyo.csv")),
                   pd.read_csv(os.path.join(lib.PROC, "validation_loyo_extra.csv"))])
    v = v[v["grid"] == 0.5]
    cols = {"modis": "MODIS input", "model": "model", "ukr_month_WAPE": "WAPE Ukraine×month",
            "ukr_year_WAPE": "WAPE Ukraine×year", "oblast_month_WAPE": "WAPE oblast×month",
            "cell_month_WAPE": "WAPE cell×month", "ukr_month_bias": "bias"}
    vt = v[list(cols)].rename(columns=cols).sort_values(["MODIS input", "model"], ascending=[False, True])
    st.dataframe(vt.style.format({c: "{:.1%}" for c in cols.values() if c.startswith(("WAPE", "bias"))}),
                 hide_index=True, use_container_width=True)
    st.caption("WAPE = Σ|pred − obs| / Σ obs. B2 = B1 × dominant land cover; M1_noFRP = GLM without FRP term. "
               "Cell = 0.5°.")
    st.markdown("**Interval calibration.** Year-bootstrap of the ratios alone gave far too narrow intervals "
                "(80% intervals covered 14% of held-out Ukraine-months). Intervals are therefore calibrated on "
                "leave-one-year-out residuals, by predicted size and AOI scale. Cross-validated coverage:")
    cv = pd.read_csv(os.path.join(lib.PROC, "validation_interval_coverage.csv"))
    st.dataframe(cv.pivot_table(index=["model", "level"], columns="interval", values="coverage")
                 .rename(columns={"80": "80% nominal", "95": "95% nominal"}).style.format("{:.0%}"),
                 use_container_width=True)
    st.markdown("**Sensitivity (fit on 2012–2025 with conflict-affected detections excluded).** "
                "22 of 24 month × day/night ratios stay inside the main fit's 95% bootstrap range "
                "(May and Sep night ratios −20%). Applying the 2012–2021 model to 2022–2025 out of sample gives "
                "+5.8% bias and 8.0% Ukraine×year WAPE. **NOAA-20 → S-NPP:** single ratio 1.011 "
                "(95% 1.000–1.023), LOYO Ukraine×month WAPE 3.7%; a month × day/night ratio did not validate better.")
    if len(AOI["cells"]):
        wk = lib.aoi_weekly(tuple(sorted(map(int, AOI["cells"]))), DOMAIN)
        yr = lib.aggregate(wk, ["year"])
        f = go.Figure()
        f.add_trace(go.Scatter(x=yr["year"], y=yr["raw_modis"], name="Raw MODIS", line=dict(color="#666666")))
        f.add_trace(go.Scatter(x=yr["year"], y=yr["raw_snpp"], name="Raw S-NPP", line=dict(color="#e66101")))
        f.add_trace(go.Scatter(x=yr["year"], y=yr["point"], name="Harmonized", line=dict(color="#31688e", width=3),
                               error_y=dict(type="data", symmetric=False, array=yr["hi95"] - yr["point"],
                                            arrayminus=yr["point"] - yr["lo95"])))
        f.update_layout(height=360, title=f"Selected AOI: {AOI['name']} ({lib.DOMAIN_LABELS[DOMAIN]})",
                        yaxis_title=UNIT, margin=dict(l=10, r=10, t=40, b=10))
        st.plotly_chart(f, use_container_width=True, config=PNG_CONFIG)


def page_methods():
    st.title("Data, methods & limitations")
    st.markdown(f"""
### What is shown
The unit everywhere is **{UNIT}**: the number of active-fire detections a VIIRS S-NPP 375 m sensor would record.
A detection is an *observation of a hot pixel*, not a fire: one fire can produce many detections, and MODIS and
VIIRS counts are never added together.

### Data
- **NASA FIRMS** active-fire archives (standard processing) and near-real-time (NRT) files for Ukraine:
  MODIS Collection 6.1 (MCD14ML / MCD14DL; Terra from Nov 2000, Aqua from Jul 2002), VIIRS 375 m S-NPP (VNP14IMG,
  2012 – 2 Nov 2026), NOAA-20 (VJ114IMG, from Apr 2018), NOAA-21 (from 2023).
- Filters: `type = 0` (presumed vegetation fire; static industrial sources removed — NRT files have no `type`,
  so a static-source mask built from standard-processing type-2 detections is applied), MODIS confidence ≥ 30,
  VIIRS confidence nominal/high.
- **ESA WorldCover 2021** (v200) land cover; **geoBoundaries** ADM0/ADM1 (Ukraine);
  Chornobyl Exclusion Zone boundary from **OpenStreetMap** (relation 3311547, © OSM contributors, ODbL).

### Conflict-affected oblasts
{lib.CONFLICT_NOTE}

### Limitations
- We count **satellite detections, not officially reported fires**. Comparison with official fire statistics
  (e.g. the State Emergency Service of Ukraine) is future work.
- No detection ≠ no fire: clouds, smoke, overpass timing and detection limits cause misses.
- Harmonized values before 2012 are **model estimates** with uncertainty (intervals shown). B1 is zero where
  MODIS saw nothing, so the *spatial* pattern of small fires that only VIIRS can see is missing before 2012
  (totals are calibrated; cell-level error ≈ 35%).
- Nov 2000 – Jul 2002 relies on Terra only (lower confidence).
- Land-cover classes come from a ~70 m overview of WorldCover; small patches are smoothed.
- All statistics are descriptive; associations are not causation; no ignition forecasts are made.
- Most recent weeks come from NRT processing, which can differ slightly from the later standard archive.
- **Anomalies and critical periods use point estimates.** Model uncertainty of the pre-2012 values is not
  propagated into percentiles or critical-week probabilities (their intervals reflect only year-to-year
  variation). Errors of the MODIS-based estimates are correlated within a year, so a proper propagation would
  need simulated whole-year error paths; this is future work. For small AOIs the default baseline is therefore
  2012–2025 (observed VIIRS only).
- **Weekly resolution:** season onset/peak/end are resolved to ±7 days; week 52 absorbs the last 8–9 days of the
  year. With ≤25 baseline years a percentile step is ~4%; weeks with zero activity in an AOI where zero is
  typical rank near the 40th–50th percentile (neutral), by design. The spring / summer-autumn split is fixed per
  AOI from its climatology, so an unusually late spring can be partly counted as summer-autumn.
- **Critical periods** use a 75th-percentile threshold pooled over all weeks of the baseline years (a per-week
  threshold would be exceeded in ~25% of years by construction); for sparse AOIs the threshold can be 0
  ("any activity").
- **Long-term trends** are computed on a constant area. The significant decline in Ukraine without
  conflict-affected oblasts (2001–2025) is not significant before the invasion (2001–2021), i.e. it depends on
  the 2022–2025 years; we do not attribute it to any cause.
""")


AOI, DOMAIN, D0, D1, BASE = sidebar()
pg = st.navigation([st.Page(page_calendar, title="Calendar", icon="📅", default=True),
                    st.Page(page_anomalies, title="Anomalies & critical periods", icon="⚠️"),
                    st.Page(page_harmonization, title="Harmonization", icon="🛰️"),
                    st.Page(page_methods, title="Data, methods & limitations", icon="📖")])
pg.run()

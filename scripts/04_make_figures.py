#!/usr/bin/env python3
"""04_make_figures.py — Generate camera-ready publication figures under `paper/figs/`.

Reads `data/processed/*.parquet` (plus working geometry in `data/processed/_work/` when
available for the raw-vs-cleaned topology zoom) and writes 300-dpi PNG figures:
  - paper/figs/fig1_network_cleaning.png
  - paper/figs/fig2_leeds_aadt_map.png
  - paper/figs/fig3_validation_scatter.png
  - paper/figs/fig4_cleaning_and_ablation.png
  - paper/figs/fig5_spatial_residuals_and_classes.png
"""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.2,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    }
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
PROC_DIR = DATA_DIR / "processed"
WORK_DIR = PROC_DIR / "_work"
FIG_DIR = REPO_ROOT / "paper" / "figs"


def _gdf_to_segments(gdf: gpd.GeoDataFrame) -> list[np.ndarray]:
    segs = []
    for geom in gdf.geometry.values:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "LineString":
            segs.append(np.asarray(geom.coords)[:, :2])
        elif geom.geom_type == "MultiLineString":
            for part in geom.geoms:
                segs.append(np.asarray(part.coords)[:, :2])
    return segs


def make_fig1_network_cleaning() -> None:
    """Fig 1: Armley Gyratory topology zoom (N1 vs N2 vs N4) + dft:95072 false-match + netclean waterfall."""
    step_df = pd.read_parquet(PROC_DIR / "netclean_step_log.parquet")
    dft_sites = pd.read_parquet(PROC_DIR / "dft_counts_leeds.parquet")

    n1 = gpd.read_parquet(WORK_DIR / "N1_raw_motor.parquet")
    n2 = gpd.read_parquet(WORK_DIR / "N2_minimal_python.parquet")
    n4 = gpd.read_parquet(WORK_DIR / "N4_csis_strategic.parquet")

    fig = plt.figure(figsize=(13.5, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.08, 0.95])

    # Top row: Armley Gyratory & A58(M) Inner Ring Road zoom (1.1 km x 1.1 km)
    xmin, xmax = 428050, 429250
    ymin, ymax = 432850, 434050
    bbox = shapely.box(xmin, ymin, xmax, ymax)

    rc_colors = {
        "motorway": "#b30000",
        "trunk": "#e34a33",
        "primary": "#fc8d59",
        "secondary": "#2b8cbe",
        "tertiary": "#7bccc4",
        "minor": "#888888",
        "service": "#cccccc",
        "other": "#aaaaaa",
    }

    for col_idx, (title, gdf) in enumerate(
        [
            ("(a) N1 Raw OSM motor (aadtbench)", n1),
            ("(b) N2 Minimal clean (Python LCC + weld)", n2),
            ("(c) N4 CSIS strategic (netclean + netweld)", n4),
        ]
    ):
        ax = fig.add_subplot(gs[0, col_idx])
        sub = gdf[gdf.intersects(bbox)].copy()
        for rc_val, grp in sub.groupby("road_class"):
            lw = 2.0 if rc_val in {"motorway", "trunk", "primary"} else (1.3 if rc_val in {"secondary", "tertiary"} else 0.7)
            lc = LineCollection(_gdf_to_segments(grp), colors=rc_colors.get(rc_val, "#888888"), linewidths=lw, alpha=0.9, zorder=2)
            ax.add_collection(lc)

        # Plot graph endpoint vertices inside zoom box
        pu = shapely.get_point(sub.geometry.values, 0)
        pv = shapely.get_point(sub.geometry.values, -1)
        pts = np.unique(
            np.round(np.c_[np.r_[shapely.get_x(pu), shapely.get_x(pv)], np.r_[shapely.get_y(pu), shapely.get_y(pv)]], 1),
            axis=0,
        )
        in_b = (pts[:, 0] >= xmin) & (pts[:, 0] <= xmax) & (pts[:, 1] >= ymin) & (pts[:, 1] <= ymax)
        pts = pts[in_b]
        ax.scatter(pts[:, 0], pts[:, 1], s=5.5, c="#111111", alpha=0.65, zorder=3, linewidths=0)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title)
        ax.text(
            0.03, 0.96,
            f"Zoom links: {len(sub):,}\nVertices: {len(pts):,}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.2,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc", alpha=0.92),
        )
        if col_idx == 0:
            # 250m scale bar in bottom-left, well clear of top-left stats box
            ax.plot([xmin + 60, xmin + 310], [ymin + 70, ymin + 70], color="#111111", lw=2.5, zorder=5)
            ax.text(xmin + 185, ymin + 95, "250 m", ha="center", va="bottom", fontsize=7.8, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.82))

    # Legend on top-right panel
    leg_items = [
        Line2D([0], [0], color="#b30000", lw=2.2, label="Motorway / Trunk"),
        Line2D([0], [0], color="#fc8d59", lw=2.0, label="Primary (A-road)"),
        Line2D([0], [0], color="#2b8cbe", lw=1.5, label="Secondary / Tertiary"),
        Line2D([0], [0], color="#888888", lw=1.0, label="Minor / Service"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#111111", markersize=4.5, label="Graph vertex"),
    ]
    fig.axes[2].legend(handles=leg_items, loc="upper right", framealpha=0.92, fontsize=7.5)

    # Bottom-left (panel d): Stourton A61 (dft:95072) industrial-yard false-match eliminated by netclean
    ax_d = fig.add_subplot(gs[1, 0])
    s95 = dft_sites[dft_sites["site_id"] == "dft:95072"].iloc[0]
    sx, sy = float(s95["x"]), float(s95["y"])
    b_d = shapely.box(sx - 280, sy - 220, sx + 220, sy + 280)
    sub_n1 = n1[n1.intersects(b_d)]
    sub_n4 = n4[n4.intersects(b_d)]
    ax_d.add_collection(LineCollection(_gdf_to_segments(sub_n1), colors="#cccccc", linewidths=1.2, linestyles="--", label="Pruned raw service stubs (N1)"))
    for rc_v, grp in sub_n4.groupby("road_class"):
        col = "#e34a33" if rc_v in {"primary", "trunk", "motorway"} else "#4a6fe3"
        lw = 2.6 if rc_v in {"primary", "trunk", "motorway"} else 1.4
        ax_d.add_collection(LineCollection(_gdf_to_segments(grp), colors=col, linewidths=lw))
    ax_d.scatter([sx], [sy], s=85, c="#d95f02", edgecolors="black", zorder=6, marker="*", label="dft:95072 ('A61', offset 147 m)")
    circle = plt.Circle((sx, sy), 30.0, fill=False, ec="#d95f02", ls=":", lw=1.5, zorder=5, label="30 m snap radius")
    ax_d.add_patch(circle)
    ax_d.set_xlim(sx - 280, sx + 220)
    ax_d.set_ylim(sy - 220, sy + 280)
    ax_d.set_aspect("equal")
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    ax_d.set_title("(d) False-snap elimination (site dft:95072, Stourton)")
    ax_d.legend(loc="lower left", framealpha=0.92, fontsize=7.4)

    # Bottom-middle & right (panel e): 14-step netclean topological reduction waterfall
    ax_e = fig.add_subplot(gs[1, 1:])
    steps = step_df["step"].tolist()
    edges_k = step_df["edges"].to_numpy(float) / 1e3
    nodes_k = step_df["nodes"].to_numpy(float) / 1e3
    comps = step_df["components"].to_numpy(int)
    x = np.arange(len(steps))
    w = 0.38
    ax_e.bar(x - w / 2, edges_k, width=w, color="#2b5c8f", label="Graph edges (thousands)")
    ax_e.bar(x + w / 2, nodes_k, width=w, color="#73a2c6", label="Graph nodes (thousands)")
    ax_e.set_xticks(x)
    ax_e.set_xticklabels(steps, rotation=32, ha="right", fontsize=8.0)
    ax_e.set_ylabel("Count (thousands)")
    ax_e.set_title("(e) 14-stage Rust netclean pipeline progression (full 15 km cordon)")
    ax_e.grid(axis="y", ls=":", alpha=0.5)

    ax_e2 = ax_e.twinx()
    ax_e2.plot(x, comps, color="#d95f02", marker="o", lw=1.8, markersize=4.5, label="Connected components")
    ax_e2.set_ylabel("Weakly connected components", color="#d95f02")
    ax_e2.tick_params(axis="y", labelcolor="#d95f02")

    h1, l1 = ax_e.get_legend_handles_labels()
    h2, l2 = ax_e2.get_legend_handles_labels()
    ax_e.legend(h1 + h2, l1 + l2, loc="upper right", framealpha=0.92)

    fig.suptitle(
        "Topological Network Cleaning in Leeds: Raw OSM vs Minimal Python vs CSIS Rust Toolchain",
        fontsize=12, fontweight="bold",
    )
    fig.text(
        0.5, -0.022,
        "Basemap data © OpenStreetMap contributors (ODbL 1.0); count site metadata © Department for Transport (OGL v3.0).",
        ha="center", fontsize=7.8, color="#444444",
    )
    fig.savefig(FIG_DIR / "fig1_network_cleaning.png", bbox_inches="tight")
    plt.close(fig)


def make_fig2_leeds_aadt_map() -> None:
    """Fig 2: Central motor traffic AADT across Leeds LAD + central core peak V/C ratio zoom."""
    lad_gdf = gpd.read_file(DATA_DIR / "leeds-lad.geojson").to_crs(27700)
    links = gpd.read_parquet(PROC_DIR / "leeds_lad_aadt.parquet")
    preds = pd.read_parquet(PROC_DIR / "site_predictions.parquet")
    p_cal = preds[(preds["network_id"] == "N4_csis_strategic") & (preds["assignment_id"] == "ue_improved_cal")].copy()

    # For dual-carriageway visual comparability with 2-way DfT counts, plot directional links at 2x or note link AADT:
    # Actually link AADT on each carriageway is the physical link flow; let's sort ascending so major corridors draw on top!
    links_sorted = links.sort_values("aadt_motor", ascending=True).reset_index(drop=True)

    bins = [0, 2000, 7500, 15000, 30000, 60000, 200000]
    labels = ["< 2k", "2k – 7.5k", "7.5k – 15k", "15k – 30k", "30k – 60k", "> 60k"]
    colors = ["#d9d9d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#67000d"]
    widths = [0.28, 0.65, 1.15, 1.80, 2.55, 3.40]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 6.8), gridspec_kw={"width_ratios": [1.18, 1.0]}, constrained_layout=True)

    # Panel (a): Full Leeds LAD central motor AADT
    lad_gdf.boundary.plot(ax=ax1, color="#222222", lw=1.2, ls="--", zorder=5)
    aadt_vals = links_sorted["aadt_motor"].to_numpy(float)
    cat_idx = np.digitize(aadt_vals, bins[1:-1])
    for bi in range(len(labels)):
        sub = links_sorted.iloc[cat_idx == bi]
        if len(sub) == 0:
            continue
        lc = LineCollection(_gdf_to_segments(sub), colors=colors[bi], linewidths=widths[bi], alpha=0.92, zorder=2 + bi * 0.1)
        ax1.add_collection(lc)

    # Overlay DfT count sites coloured by observed AADF (halved when on dual carriageway? No, plot observed 2-way with matching bins!)
    obs_cat = np.digitize(p_cal["obs_motor"].to_numpy(float), bins[1:-1])
    for bi in range(len(labels)):
        sp = p_cal.iloc[obs_cat == bi]
        if len(sp) == 0:
            continue
        ax1.scatter(sp["x"], sp["y"], s=14 + bi * 5, c=colors[bi], edgecolors="#111111", linewidths=0.55, zorder=6)

    minx, miny, maxx, maxy = lad_gdf.total_bounds
    pad = 1200
    ax1.set_xlim(minx - pad, maxx + pad)
    ax1.set_ylim(miny - pad, maxy + pad)
    ax1.set_aspect("equal")
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlabel("")
    ax1.set_ylabel("")
    ax1.set_title("(a) Modelled Motor Traffic AADT & DfT Count Sites across Leeds LAD")

    # 5 km scale bar
    ax1.plot([minx + 800, minx + 5800], [miny + 800, miny + 800], color="#111111", lw=2.8, zorder=7)
    ax1.text(minx + 3300, miny + 1250, "5 km", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    leg_a = [Line2D([0], [0], color=colors[i], lw=widths[i] + 0.5, label=f"{labels[i]} veh/day") for i in range(len(labels))]
    leg_a.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#fb6a4a", markeredgecolor="#111111", markersize=6, label="DfT count site (2021–25)"))
    leg_a.append(Line2D([0], [0], color="#222222", lw=1.2, ls="--", label="Leeds LAD (E08000035)"))
    ax1.legend(handles=leg_a, loc="lower right", framealpha=0.94, title="Daily Motor AADT", title_fontsize=8.5)

    # Panel (b): Urban Core Zoom (7 km x 7 km) showing Peak-Hour V/C Ratio from Multi-Class BFW UE
    zx0, zx1 = 426000, 433500
    zy0, zy1 = 429800, 437300
    zbox = shapely.box(zx0, zy0, zx1, zy1)
    core = links[links.intersects(zbox)].sort_values("vc_ratio_peak", ascending=True).copy()

    vc_bins = [0.0, 0.25, 0.50, 0.75, 0.95, 10.0]
    vc_labels = ["< 0.25 (Free flow)", "0.25 – 0.50 (Light)", "0.50 – 0.75 (Moderate)", "0.75 – 0.95 (Busy)", ">= 0.95 (At/over cap.)"]
    vc_cols = ["#d1e5f0", "#67a9cf", "#fdcc8a", "#ef6548", "#990000"]
    vc_w = [0.45, 1.0, 1.6, 2.3, 3.0]
    vc_idx = np.digitize(core["vc_ratio_peak"].to_numpy(float), vc_bins[1:-1])
    for bi in range(len(vc_labels)):
        sub = core.iloc[vc_idx == bi]
        if len(sub) == 0:
            continue
        ax2.add_collection(LineCollection(_gdf_to_segments(sub), colors=vc_cols[bi], linewidths=vc_w[bi], alpha=0.92))

    ax2.set_xlim(zx0, zx1)
    ax2.set_ylim(zy0, zy1)
    ax2.set_aspect("equal")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title("(b) Central Leeds & M621 Corridor: Peak-Hour V/C Ratio (UE-Improved)")

    # Annotate key corridors
    annots = [
        (429800, 431200, "M621"),
        (429600, 434150, "A58(M) IRR"),
        (428400, 433400, "Armley Gyratory"),
        (431800, 434500, "A64 York Rd"),
        (428100, 436100, "A660 Headingley"),
    ]
    for ax_x, ax_y, txt in annots:
        ax2.text(
            ax_x, ax_y, txt, fontsize=7.6, fontweight="bold", color="#111111",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.78),
        )

    ax2.plot([zx0 + 400, zx0 + 1400], [zy0 + 400, zy0 + 400], color="#111111", lw=2.6)
    ax2.text(zx0 + 900, zy0 + 550, "1 km", ha="center", va="bottom", fontsize=8.2, fontweight="bold")
    leg_b = [Line2D([0], [0], color=vc_cols[i], lw=vc_w[i] + 0.4, label=vc_labels[i]) for i in range(len(vc_labels))]
    ax2.legend(handles=leg_b, loc="lower right", framealpha=0.94, title="Peak-Hour PCE V/C", title_fontsize=8.5)

    # Draw zoom box on Panel (a)
    ax1.plot([zx0, zx1, zx1, zx0, zx0], [zy0, zy0, zy1, zy1, zy0], color="#1f78b4", lw=1.1, ls=":", zorder=7)

    fig.text(
        0.5, -0.018,
        "Contains OpenStreetMap data © OpenStreetMap contributors (ODbL 1.0); ONS LAD boundary & DfT AADF counts under Open Government Licence v3.0.",
        ha="center", fontsize=8.0, color="#444444",
    )
    fig.savefig(FIG_DIR / "fig2_leeds_aadt_map.png", bbox_inches="tight")
    plt.close(fig)


def _geh_envelope(y_grid: np.ndarray, geh_thresh: float = 5.0) -> tuple[np.ndarray, np.ndarray]:
    """Solve 2*(0.1*p - 0.1*y)^2 / (0.1*p + 0.1*y) = G^2 for upper & lower daily AADT bounds."""
    # Let h = 0.1*y, x = 0.1*p. 2(x - h)^2 = G^2 (x + h) => 2 x^2 - (4h + G^2) x + (2h^2 - G^2 h) = 0
    h = 0.1 * y_grid
    g2 = geh_thresh ** 2
    b = 4.0 * h + g2
    disc = np.maximum(b ** 2 - 8.0 * (2.0 * h ** 2 - g2 * h), 0.0)
    x_lo = np.maximum((b - np.sqrt(disc)) / 4.0, 1.0)
    x_hi = (b + np.sqrt(disc)) / 4.0
    return x_lo * 10.0, x_hi * 10.0


def make_fig3_validation_scatter() -> None:
    """Fig 3: 2x2 Observed vs Predicted AADT scatter with WebTAG GEH < 5 envelope."""
    preds = pd.read_parquet(PROC_DIR / "site_predictions.parquet")
    metrics = pd.read_parquet(PROC_DIR / "validation_metrics.parquet")
    m_lad = metrics[metrics["eval_set"] == "leeds_lad_2021_2025"]

    panels = [
        ("N0_raw_unfiltered", "aon_physical", "(a) Uncritical Baseline: Raw OSM + AON"),
        ("N1_raw_motor", "ue_bfw_physical", "(b) Benchmark Baseline: Raw Motor + UE-BFW"),
        ("N4_csis_strategic", "ue_improved_phys", "(c) Cleaned + Physical Multi-Class UE (k = 1)"),
        ("N4_csis_strategic", "ue_improved_cal", "(d) Cleaned + CSIS Out-of-Fold Calibrated"),
    ]
    cat_palette = {
        "TM": ("#b30000", "Trunk Motorway (TM)"),
        "PM": ("#e34a33", "Principal Motorway (PM)"),
        "TA": ("#fe9929", "Trunk A-road (TA)"),
        "PA": ("#2b8cbe", "Principal A-road (PA)"),
        "MB": ("#41ab5d", "Minor B-road (MB)"),
        "MCU": ("#756bb1", "Minor C/Unclassified (MCU)"),
    }

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.8), constrained_layout=True)
    y_grid = np.logspace(2.5, 5.35, 300)
    lo5, hi5 = _geh_envelope(y_grid, 5.0)

    for ax, (nid, aid, title) in zip(axes.ravel(), panels):
        sub = preds[(preds["network_id"] == nid) & (preds["assignment_id"] == aid)].copy()
        mrow = m_lad[(m_lad["network_id"] == nid) & (m_lad["assignment_id"] == aid)].iloc[0]

        ax.fill_between(y_grid, lo5, hi5, color="#2ca25f", alpha=0.14, label="WebTAG GEH < 5 band")
        ax.plot(y_grid, y_grid, color="#222222", lw=1.2, ls="--", label="1:1 agreement")

        obs = np.maximum(sub["obs_motor"].to_numpy(float), 250.0)
        prd = np.maximum(sub["pred_motor"].to_numpy(float), 250.0)
        for cat_code, (col, lbl) in cat_palette.items():
            mask = (sub["road_category"] == cat_code).to_numpy(bool)
            if mask.any():
                ax.scatter(obs[mask], prd[mask], s=26, c=col, alpha=0.82, edgecolors="white", linewidths=0.35, label=lbl)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(300, 220000)
        ax.set_ylim(200, 250000)
        ax.grid(True, which="major", ls=":", alpha=0.5)
        ax.set_title(title)
        ax.set_xlabel("Observed DfT Motor AADF (veh/day, 2021–2025)")
        ax.set_ylabel("Modelled Motor AADT (veh/day)")

        n_zero = int((sub["pred_motor"] < 260.0).sum())
        stat_txt = (
            f"Matched n = {int(mrow['n_sites_matched'])} (≤250 clamp: {n_zero})\n"
            f"Pearson r = {mrow['pearson_r']:.3f} | Spearman ρ = {mrow['spearman_rho']:.3f}\n"
            f"MAE = {mrow['mae']:,.0f} | MedAPE = {mrow['median_ape_pct']:.1f}%\n"
            f"GEH < 5 = {mrow['geh_lt_5_pct']:.1f}% | GEH < 10 = {mrow['geh_lt_10_pct']:.1f}%"
        )
        ax.text(
            0.03, 0.96, stat_txt, transform=ax.transAxes, va="top", ha="left", fontsize=8.1,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.93),
        )

    handles, labels_leg = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels_leg, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.055), frameon=True)
    fig.suptitle("Out-of-Sample Validation Against Leeds LAD DfT Motor Counts (5-Fold 3 km Spatial CV)", fontsize=12, fontweight="bold")
    fig.savefig(FIG_DIR / "fig3_validation_scatter.png", bbox_inches="tight")
    plt.close(fig)


def make_fig4_cleaning_and_ablation() -> None:
    """Fig 4: Network cleaning sensitivity across assignment methods + 7-step ablation & MVKM check."""
    metrics = pd.read_parquet(PROC_DIR / "validation_metrics.parquet")
    vkm_ref = pd.read_parquet(PROC_DIR / "dft_vkm_reference.parquet")
    m_lad = metrics[metrics["eval_set"] == "leeds_lad_2021_2025"].copy()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.5, 4.9), constrained_layout=True)

    # Panel (a): Cleaning progression N0..N4 across assignment variants
    net_order = ["N0_raw_unfiltered", "N1_raw_motor", "N2_minimal_python", "N3_csis_netclean", "N4_csis_strategic"]
    net_short = ["N0\nRaw all", "N1\nRaw motor", "N2\nMin Python", "N3\nCSIS full", "N4\nCSIS strat."]
    series_specs = [
        ("aon_physical", "AON (physical)", "#7f7f7f", "o", "--"),
        ("ue_bfw_physical", "UE-BFW (physical baseline)", "#e66101", "s", ":"),
        ("ue_bfw_aadtbench", "UE-BFW (aadtbench modulated)", "#fdb863", "^", "-."),
        ("ue_improved_phys", "UE-Improved (physical k=1)", "#2b8cbe", "D", "-"),
        ("ue_improved_cal", "UE-Improved + CSIS OOF cal.", "#084081", "*", "-"),
    ]
    for aid, lbl, col, mk, ls in series_specs:
        xs, ys = [], []
        for idx_n, nid in enumerate(net_order):
            r = m_lad[(m_lad["network_id"] == nid) & (m_lad["assignment_id"] == aid)]
            if not r.empty:
                xs.append(idx_n)
                ys.append(float(r["pearson_r"].iloc[0]))
        ax1.plot(xs, ys, label=lbl, color=col, marker=mk, ls=ls, lw=2.0, markersize=6.5)

    ax1.set_xticks(range(len(net_order)))
    ax1.set_xticklabels(net_short)
    ax1.set_ylim(0.45, 0.94)
    ax1.set_ylabel("Linear Pearson correlation r vs DfT AADF")
    ax1.set_title("(a) Benefit of CSIS Cleaning Across Methods")
    ax1.grid(True, ls=":", alpha=0.55)
    ax1.legend(loc="lower right", fontsize=7.6, framealpha=0.94)

    # Panel (b): 7-Stage Methodological Ablation (Pearson r, Spearman rho, GEH<5%, GEH<10%)
    ablation_steps = [
        ("N0_raw_unfiltered", "aon_physical", "A0: Raw all\n+ AON"),
        ("N1_raw_motor", "aon_physical", "A1: Motor\nfilter + AON"),
        ("N1_raw_motor", "ue_bfw_physical", "A2: Raw\n+ UE-BFW"),
        ("N4_csis_strategic", "ue_bfw_physical", "A3: CSIS N4\n+ UE-BFW"),
        ("N4_csis_strategic", "ue_welded_internal", "A4: + netweld\nlanes/speed"),
        ("N4_csis_strategic", "ue_improved_phys_scaled", "A5: + SRN &\nmulti-class"),
        ("N4_csis_strategic", "ue_improved_cal", "A6: + CSIS\nOOF cal."),
    ]
    ab_rows = []
    for nid, aid, albl in ablation_steps:
        r = m_lad[(m_lad["network_id"] == nid) & (m_lad["assignment_id"] == aid)].iloc[0]
        ab_rows.append({"lbl": albl, **r.to_dict()})
    ab_df = pd.DataFrame(ab_rows)

    x_ab = np.arange(len(ab_df))
    ax2.bar(x_ab - 0.18, ab_df["geh_lt_5_pct"], width=0.35, color="#2ca25f", label="GEH < 5 (%)")
    ax2.bar(x_ab + 0.18, ab_df["geh_lt_10_pct"], width=0.35, color="#99d8c9", label="GEH < 10 (%)")
    ax2.set_xticks(x_ab)
    ax2.set_xticklabels(ab_df["lbl"], rotation=28, ha="right", fontsize=7.8)
    ax2.set_ylabel("Share of DfT sites (%)")
    ax2.set_ylim(0, 62)
    ax2.set_title("(b) Stepwise Ablation: WebTAG GEH Pass Rate")
    ax2.grid(axis="y", ls=":", alpha=0.55)
    ax2.legend(loc="upper left", fontsize=8.0)

    # Panel (c): Network-level Leeds LAD Daily Motor-Vehicle-Km (MVKM) vs DfT Reference
    maj_ref = float(vkm_ref.loc[vkm_ref["road_category"] == "ALL_MAJOR", "mvkm_day_postcovid"].iloc[0])
    all_ref = float(vkm_ref.loc[vkm_ref["road_category"] == "IMPLIED_ALL_ROADS_DFT_SHARE", "mvkm_day_postcovid"].iloc[0])

    minor_mvkm = ab_df["lad_mvkm_day_all"] - ab_df["lad_mvkm_day_major"]
    ax3.bar(x_ab, ab_df["lad_mvkm_day_major"], width=0.55, color="#d95f02", label="Modelled Major Roads (M/A)")
    ax3.bar(x_ab, np.maximum(minor_mvkm, 0.0), bottom=ab_df["lad_mvkm_day_major"], width=0.55, color="#7570b3", alpha=0.82, label="Modelled Minor Roads (B/C/U)")
    ax3.axhline(maj_ref, color="#b30000", ls="--", lw=1.6, label=f"DfT Major Road ref. ({maj_ref:.1f}M)")
    ax3.axhline(all_ref, color="#222222", ls=":", lw=1.6, label=f"DfT Implied All-Roads ({all_ref:.1f}M)")
    ax3.set_xticks(x_ab)
    ax3.set_xticklabels(ab_df["lbl"], rotation=28, ha="right", fontsize=7.8)
    ax3.set_ylabel("Leeds LAD Million Motor-Veh-Km / Day")
    ax3.set_title("(c) Network Total Traffic Volume vs DfT")
    ax3.grid(axis="y", ls=":", alpha=0.55)
    ax3.legend(loc="upper right", fontsize=7.5, framealpha=0.94)

    fig.savefig(FIG_DIR / "fig4_cleaning_and_ablation.png", bbox_inches="tight")
    plt.close(fig)


def make_fig5_spatial_residuals_and_classes() -> None:
    """Fig 5: Spatial distribution of WebTAG GEH across Leeds LAD + road-category mean AADT comparison."""
    lad_gdf = gpd.read_file(DATA_DIR / "leeds-lad.geojson").to_crs(27700)
    links = gpd.read_parquet(PROC_DIR / "leeds_lad_aadt.parquet")
    preds = pd.read_parquet(PROC_DIR / "site_predictions.parquet")
    cat_df = pd.read_parquet(PROC_DIR / "category_breakdown.parquet")

    p_cal = preds[(preds["network_id"] == "N4_csis_strategic") & (preds["assignment_id"] == "ue_improved_cal")].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.8), gridspec_kw={"width_ratios": [1.1, 1.0]}, constrained_layout=True)

    # Panel (a): Spatial GEH map on Leeds LAD
    lad_gdf.boundary.plot(ax=ax1, color="#222222", lw=1.1, ls="--", zorder=2)
    maj_links = links[links["road_class"].isin(["motorway", "trunk", "primary", "secondary"])]
    ax1.add_collection(LineCollection(_gdf_to_segments(maj_links), colors="#cccccc", linewidths=0.85, zorder=1))

    g_ok = p_cal[p_cal["geh"] < 5.0]
    g_mid = p_cal[(p_cal["geh"] >= 5.0) & (p_cal["geh"] < 10.0)]
    g_hi = p_cal[p_cal["geh"] >= 10.0]

    ax1.scatter(g_ok["x"], g_ok["y"], s=38, c="#2ca25f", marker="o", edgecolors="#111111", linewidths=0.5, zorder=5, label=f"GEH < 5 (n={len(g_ok)}, {len(g_ok)/len(p_cal)*100:.1f}%)")
    ax1.scatter(g_mid["x"], g_mid["y"], s=42, c="#fe9929", marker="^", edgecolors="#111111", linewidths=0.5, zorder=4, label=f"5 ≤ GEH < 10 (n={len(g_mid)}, {len(g_mid)/len(p_cal)*100:.1f}%)")
    ax1.scatter(g_hi["x"], g_hi["y"], s=36, c="#de2d26", marker="s", alpha=0.82, edgecolors="#111111", linewidths=0.45, zorder=3, label=f"GEH ≥ 10 (n={len(g_hi)}, {len(g_hi)/len(p_cal)*100:.1f}%)")

    minx, miny, maxx, maxy = lad_gdf.total_bounds
    ax1.set_xlim(minx - 1000, maxx + 1000)
    ax1.set_ylim(miny - 1000, maxy + 1000)
    ax1.set_aspect("equal")
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlabel("")
    ax1.set_ylabel("")
    ax1.set_title("(a) Spatial Distribution of Out-of-Fold WebTAG GEH (N4 + Calibrated UE)")
    ax1.legend(loc="lower right", framealpha=0.94, fontsize=8.2)

    # Panel (b): Mean AADT by DfT Road Category
    cats = ["TM", "PM", "PA", "TA", "MB", "MCU"]
    cat_labels = ["TM\n(Trunk M-way)", "PM\n(Princ. M-way)", "PA\n(Princ. A-road)", "TA\n(Trunk A-road)", "MB\n(Minor B-road)", "MCU\n(C / Unclass.)"]

    c_aon = cat_df[(cat_df["network_id"] == "N1_raw_motor") & (cat_df["assignment_id"] == "aon_physical")].set_index("road_category")
    c_ue0 = cat_df[(cat_df["network_id"] == "N1_raw_motor") & (cat_df["assignment_id"] == "ue_bfw_physical")].set_index("road_category")
    c_phys = cat_df[(cat_df["network_id"] == "N4_csis_strategic") & (cat_df["assignment_id"] == "ue_improved_phys")].set_index("road_category")
    c_cal = cat_df[(cat_df["network_id"] == "N4_csis_strategic") & (cat_df["assignment_id"] == "ue_improved_cal")].set_index("road_category")

    xc = np.arange(len(cats))
    bw = 0.16
    ax2.bar(xc - 2 * bw, [c_cal.loc[c, "obs_mean"] / 1e3 for c in cats], width=bw, color="#222222", label="Observed DfT AADF")
    ax2.bar(xc - 1 * bw, [c_aon.loc[c, "pred_mean"] / 1e3 for c in cats], width=bw, color="#bdbdbd", label="N1 + AON")
    ax2.bar(xc + 0 * bw, [c_ue0.loc[c, "pred_mean"] / 1e3 for c in cats], width=bw, color="#fdb863", label="N1 + UE-BFW")
    ax2.bar(xc + 1 * bw, [c_phys.loc[c, "pred_mean"] / 1e3 for c in cats], width=bw, color="#2b8cbe", label="N4 + UE-Improved (phys)")
    ax2.bar(xc + 2 * bw, [c_cal.loc[c, "pred_mean"] / 1e3 for c in cats], width=bw, color="#084081", label="N4 + UE-Improved (OOF cal)")

    ax2.set_xticks(xc)
    ax2.set_xticklabels(cat_labels, fontsize=8.0)
    ax2.set_ylabel("Mean Motor AADT (thousand veh/day)")
    ax2.set_title("(b) Hierarchy Realism Across DfT Road Categories")
    ax2.grid(axis="y", ls=":", alpha=0.55)
    ax2.legend(loc="upper right", fontsize=8.0, framealpha=0.94)

    fig.text(
        0.5, -0.022,
        "Contains OpenStreetMap data © OpenStreetMap contributors (ODbL 1.0); ONS & DfT data under Open Government Licence v3.0.",
        ha="center", fontsize=7.8, color="#444444",
    )
    fig.savefig(FIG_DIR / "fig5_spatial_residuals_and_classes.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("Generating Fig 1 (Network cleaning zoom + false-snap + waterfall)...")
    make_fig1_network_cleaning()
    print("Generating Fig 2 (Leeds LAD central AADT map + core V/C zoom)...")
    make_fig2_leeds_aadt_map()
    print("Generating Fig 3 (Observed vs predicted AADT scatter + GEH < 5 envelope)...")
    make_fig3_validation_scatter()
    print("Generating Fig 4 (Network cleaning sensitivity + ablation + MVKM check)...")
    make_fig4_cleaning_and_ablation()
    print("Generating Fig 5 (Spatial GEH map + road category breakdown)...")
    make_fig5_spatial_residuals_and_classes()
    print("All 5 figures saved to paper/figs/.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""02_clean_networks.py — Build and compare OSM network cleaning variants for Leeds.

Constructs five progressive network variants over the 15 km Leeds cordon (`leeds-v1`)
and clips/summarises each inside the ONS Leeds LAD boundary (`E08000035`):

  1. `N0_raw_unfiltered`: Raw OSM (`aadtbench` pyosmium graph, 230,330 links) with no
     modal access filtering (`allow_car=True` forced on all highway ways) — uncritical baseline.
  2. `N1_raw_motor`: Raw OSM motor-accessible links (`allow_car == True`, 143,151 links)
     from `leeds-v1/network.parquet` with zero topological cleaning.
  3. `N2_minimal_python`: Minimal topological cleaning in Python (`geopandas` + `scipy` +
     `shapely.line_merge`): motor-accessible filter, largest weakly connected component,
     iterative dead-end trimming (service & <250 m minor stubs inside 14.2 km cordon),
     and homogeneous degree-2 node contraction.
  4. `N3_csis_netclean`: Full open Rust toolchain (`osmget` -> `osmactive` GB rules ->
     `netfit network` [`netclean` 14-step consolidation, junction clustering, parallel-way
     merging, despining, and barrier splitting] -> `osmactive continuity` -> `netweld join`
     transferring parsed directional lanes and posted `maxspeed`).
  5. `N4_csis_strategic`: Strategic motor network derived from `N3_csis_netclean` by
     pruning private/parking `service`/`track` ways, severing graph connectivity across
     `osmactive`-classified motor-blocking modal filters (`allow_car == False` barriers),
     and extracting the single weakly connected motor component.

Committed outputs (<1 MB):
  - data/processed/network_variants_summary.parquet
  - data/processed/netclean_step_log.parquet

Working assignment networks (gitignored in data/processed/_work/):
  - data/processed/_work/{N0_raw_unfiltered,N1_raw_motor,N2_minimal_python,N3_csis_netclean,N4_csis_strategic}.parquet
"""
from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Transformer
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
OUT_DIR = DATA_DIR / "processed"
WORK_DIR = OUT_DIR / "_work"

CASE_DIR = Path("/home/rbn/github/robinlovelace/aadtbench/cases/leeds-v1")
OSMGET_BIN = Path("/home/rbn/github/robinlovelace/osmget/target/release/osmget")
OSMACTIVE_BIN = Path("/home/rbn/github/atf-lab/osmactive/target/release/osmactive")
OSMACTIVE_RULES = Path("/home/rbn/github/atf-lab/osmactive/rules/gb.toml")
NETFIT_BIN = Path("/home/rbn/github/atf-lab/netfit/target/release/netfit")
NETWELD_BIN = Path("/home/rbn/github/atfutures/netweld/target/release/netweld")

CENTRE_LON, CENTRE_LAT = -1.5437941, 53.7974185
RADIUS_M = 15000

HIGHWAY_TO_CLASS = {
    "motorway": "motorway", "motorway_link": "motorway",
    "trunk": "trunk", "trunk_link": "trunk",
    "primary": "primary", "primary_link": "primary",
    "secondary": "secondary", "secondary_link": "secondary",
    "tertiary": "tertiary", "tertiary_link": "tertiary",
    "unclassified": "minor", "residential": "minor", "living_street": "minor",
    "road": "minor", "connector": "minor", "service": "service",
    "cycleway": "cycleway", "busway": "other",
    "footway": "footway", "pedestrian": "footway", "steps": "footway",
    "path": "footway", "bridleway": "footway", "track": "service",
}

DEFAULT_DIR_LANES = {
    "motorway": 3.0, "motorway_link": 1.3,
    "trunk": 2.0, "trunk_link": 1.2,
    "primary": 1.25, "primary_link": 1.0,
    "secondary": 1.10, "secondary_link": 1.0,
    "tertiary": 1.0, "tertiary_link": 1.0,
    "unclassified": 0.9, "residential": 0.8, "living_street": 0.5,
    "connector": 0.9, "service": 0.4, "track": 0.3,
}


def git_sha(repo_path: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo_path), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return "unknown"


def parse_maxspeed_kmh(val: object) -> float:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return np.nan
    s = str(val).strip().lower()
    if "nsl_dual" in s:
        return 112.65
    if "nsl_single" in s or "national" in s:
        return 96.56
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    if not m:
        return np.nan
    num = float(m.group(1))
    if "mph" in s:
        return num * 1.609344
    return num


def parse_dir_lanes(ways_df: pd.DataFrame) -> np.ndarray:
    """Compute per-direction lane count from OSM `lanes`, `lanes:forward`, `lanes:backward`, and `dir_car`."""
    hw = ways_df["highway"].fillna("unclassified").astype(str)
    defaults = hw.map(DEFAULT_DIR_LANES).fillna(1.0).to_numpy(float)
    dcar = ways_df["dir_car"].fillna("both").astype(str).to_numpy()
    is_oneway = ~np.isin(dcar, ["both", "0"])

    nan_s = pd.Series(np.nan, index=ways_df.index)
    lanes_raw = pd.to_numeric(ways_df["lanes"] if "lanes" in ways_df.columns else nan_s, errors="coerce").to_numpy(float)
    lf_col = "lanes_forward" if "lanes_forward" in ways_df.columns else ("lanes:forward" if "lanes:forward" in ways_df.columns else None)
    lb_col = "lanes_backward" if "lanes_backward" in ways_df.columns else ("lanes:backward" if "lanes:backward" in ways_df.columns else None)
    lf = pd.to_numeric(ways_df[lf_col] if lf_col else nan_s, errors="coerce").to_numpy(float)
    lb = pd.to_numeric(ways_df[lb_col] if lb_col else nan_s, errors="coerce").to_numpy(float)

    out = defaults.copy()
    has_tot = np.isfinite(lanes_raw) & (lanes_raw > 0)
    out = np.where(has_tot, np.where(is_oneway, lanes_raw, lanes_raw / 2.0), out)
    has_dir = np.isfinite(lf) & np.isfinite(lb) & (lf + lb > 0)
    out = np.where(has_dir & (~is_oneway), 0.5 * (lf + lb), out)
    return np.clip(out, 0.4, 6.0)


def assign_endpoint_nodes(gdf: gpd.GeoDataFrame, blocked_xy: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Assign 1-based endpoint node IDs from 0.1 m-snapped coordinates, optionally severing blocked barrier nodes."""
    pts_u = shapely.get_point(gdf.geometry.values, 0)
    pts_v = shapely.get_point(gdf.geometry.values, -1)
    coords = np.c_[
        np.round(np.r_[shapely.get_x(pts_u), shapely.get_x(pts_v)], 1),
        np.round(np.r_[shapely.get_y(pts_u), shapely.get_y(pts_v)], 1),
    ]
    uniq_xy, inv = np.unique(coords, axis=0, return_inverse=True)
    n = len(gdf)
    u_ids = inv[:n] + 1
    v_ids = inv[n:] + 1

    if blocked_xy is not None and len(blocked_xy) > 0:
        d_bar, _ = cKDTree(blocked_xy).query(uniq_xy)
        blocked_nodes = set(np.where(d_bar <= 1.5)[0] + 1)
        next_id = len(uniq_xy) + 1
        for i in range(n):
            if int(u_ids[i]) in blocked_nodes:
                u_ids[i] = next_id
                next_id += 1
            if int(v_ids[i]) in blocked_nodes:
                v_ids[i] = next_id
                next_id += 1
    return u_ids.astype(np.int64), v_ids.astype(np.int64)


def build_minimal_clean(raw_motor: gpd.GeoDataFrame, cx: float, cy: float) -> gpd.GeoDataFrame:
    """N2: Pure-Python motor LCC + cordon-safe dangling-link trimming + homogeneous degree-2 welding."""
    df = raw_motor.copy().reset_index(drop=True)

    # 1. Largest weakly connected component
    nodes = np.unique(df[["u", "v"]].values)
    idx = pd.Series(np.arange(len(nodes)), index=nodes)
    g = coo_matrix((np.ones(len(df)), (idx[df["u"]].values, idx[df["v"]].values)), shape=(len(nodes), len(nodes)))
    _, lab = connected_components(g, directed=True, connection="weak")
    lcc_nodes = set(nodes[lab == np.bincount(lab).argmax()])
    df = df[df["u"].isin(lcc_nodes) & df["v"].isin(lcc_nodes)].copy().reset_index(drop=True)

    # Node coordinates for cordon check (protect perimeter arterial portals at r >= 14.2 km)
    pts_u = shapely.get_point(df.geometry.values, 0)
    pts_v = shapely.get_point(df.geometry.values, -1)
    n_xy = pd.concat(
        [
            pd.DataFrame({"node": df["u"].values, "x": shapely.get_x(pts_u), "y": shapely.get_y(pts_u)}),
            pd.DataFrame({"node": df["v"].values, "x": shapely.get_x(pts_v), "y": shapely.get_y(pts_v)}),
        ]
    ).drop_duplicates("node").set_index("node")
    n_xy["r"] = np.hypot(n_xy["x"].values - cx, n_xy["y"].values - cy)
    perimeter_nodes = set(n_xy.index[n_xy["r"] >= 14200.0])

    # 2. Iterative dangling-link trimming (dead-end leaves: service of any length, minor/other < 250 m)
    for _ in range(12):
        deg = pd.Series(np.r_[df["u"].values, df["v"].values]).value_counts()
        u_deg = df["u"].map(deg).values
        v_deg = df["v"].map(deg).values
        u_leaf = (u_deg == 1) & (~df["u"].isin(perimeter_nodes).values)
        v_leaf = (v_deg == 1) & (~df["v"].isin(perimeter_nodes).values)
        is_leaf = u_leaf | v_leaf
        trimmable = (df["road_class"] == "service").values | (
            df["road_class"].isin(["minor", "other"]).values & (df["length_m"].values < 250.0)
        )
        drop = is_leaf & trimmable
        if not drop.any():
            break
        df = df.loc[~drop].reset_index(drop=True)

    # 3. Homogeneous degree-2 node contraction via shapely.line_merge per attribute group
    # Identify nodes with incident link count == 2 whose two incident links belong to the same attribute key
    df["_key"] = (
        df["road_class"].fillna("") + "|"
        + df["dir_car"].fillna("both").astype(str) + "|"
        + df["allow_heavy"].fillna(True).astype(str) + "|"
        + df["name"].fillna("") + "|"
        + df["ref"].fillna("")
    )
    inc_count = pd.Series(np.r_[df["u"].values, df["v"].values]).value_counts()
    deg2_nodes = set(inc_count.index[inc_count == 2]) - perimeter_nodes

    # Junction nodes are any non-deg2 node OR any deg2 node where incident keys differ
    end_keys = pd.concat(
        [
            pd.DataFrame({"node": df["u"].values, "k": df["_key"].values}),
            pd.DataFrame({"node": df["v"].values, "k": df["_key"].values}),
        ]
    )
    n_k = end_keys.groupby("node")["k"].nunique()
    weldable = deg2_nodes & set(n_k.index[n_k == 1])

    # Split into links not touching weldable nodes vs candidate chains
    touches_weld = df["u"].isin(weldable) | df["v"].isin(weldable)
    keep_direct = df.loc[~touches_weld].drop(columns=["_key"]).copy()
    cand = df.loc[touches_weld].copy()

    merged_rows = []
    for key_val, grp in cand.groupby("_key", sort=False):
        rc, dcar, ahvy, nm, rf = key_val.split("|", 4)
        # Split multiline at any non-weldable junction node by slightly Perturbing? Even simpler:
        # Build connected components of links in `grp` sharing ONLY `weldable` nodes!
        m = len(grp)
        if m == 1:
            merged_rows.append(grp.drop(columns=["_key"]))
            continue
        # Map shared weldable nodes -> connect links sharing that weldable node
        node_to_link: dict[int, list[int]] = {}
        u_arr, v_arr = grp["u"].values, grp["v"].values
        for li in range(m):
            if int(u_arr[li]) in weldable:
                node_to_link.setdefault(int(u_arr[li]), []).append(li)
            if int(v_arr[li]) in weldable:
                node_to_link.setdefault(int(v_arr[li]), []).append(li)
        rows_i, cols_i = [], []
        for llist in node_to_link.values():
            if len(llist) == 2:
                rows_i.extend([llist[0], llist[1]])
                cols_i.extend([llist[1], llist[0]])
        if not rows_i:
            merged_rows.append(grp.drop(columns=["_key"]))
            continue
        cg = coo_matrix((np.ones(len(rows_i)), (rows_i, cols_i)), shape=(m, m))
        n_cc, cc_lab = connected_components(cg, directed=False)
        geoms = grp.geometry.values
        seg_ids = grp["segment_id"].values
        new_recs = []
        for cid in range(n_cc):
            idx_c = np.where(cc_lab == cid)[0]
            if len(idx_c) == 1:
                row0 = grp.iloc[idx_c[0]].drop(labels=["_key"]).to_dict()
                new_recs.append(row0)
            else:
                lm = shapely.line_merge(shapely.multilinestrings(geoms[idx_c]))
                parts = list(lm.geoms) if lm.geom_type == "MultiLineString" else [lm]
                for pi, part in enumerate(parts):
                    # Preserve direction along first matching sub-segment
                    p0 = shapely.get_point(part, 0)
                    # Check if p0 matches start of any sub-segment or end of all
                    starts = shapely.get_point(geoms[idx_c], 0)
                    if np.min(shapely.distance(starts, p0)) > 0.5:
                        part = shapely.reverse(part)
                    new_recs.append(
                        {
                            "segment_id": f"{seg_ids[idx_c[0]]}_m{pi}",
                            "osm_id": int(grp.iloc[idx_c[0]]["osm_id"]) if "osm_id" in grp.columns else 0,
                            "u": 0,
                            "v": 0,
                            "length_m": float(part.length),
                            "road_class": rc,
                            "highway": grp.iloc[idx_c[0]].get("highway", rc),
                            "name": nm if nm else None,
                            "ref": rf if rf else None,
                            "allow_walking": bool(grp.iloc[idx_c[0]].get("allow_walking", True)),
                            "allow_cycling": bool(grp.iloc[idx_c[0]].get("allow_cycling", True)),
                            "allow_car": True,
                            "allow_heavy": ahvy == "True",
                            "dir_cycling": "both",
                            "dir_car": dcar,
                            "geometry": part,
                        }
                    )
        merged_rows.append(gpd.GeoDataFrame(new_recs, geometry="geometry", crs=df.crs))

    out = gpd.GeoDataFrame(pd.concat([keep_direct] + merged_rows, ignore_index=True), geometry="geometry", crs=df.crs)
    out["length_m"] = out.geometry.length
    u_new, v_new = assign_endpoint_nodes(out)
    out["u"] = u_new
    out["v"] = v_new
    # Re-ensure single LCC after contraction
    nodes2 = np.unique(out[["u", "v"]].values)
    idx2 = pd.Series(np.arange(len(nodes2)), index=nodes2)
    g2 = coo_matrix((np.ones(len(out)), (idx2[out["u"]].values, idx2[out["v"]].values)), shape=(len(nodes2), len(nodes2)))
    _, lab2 = connected_components(g2, directed=True, connection="weak")
    lcc2 = set(nodes2[lab2 == np.bincount(lab2).argmax()])
    out = out[out["u"].isin(lcc2) & out["v"].isin(lcc2)].reset_index(drop=True)
    return out


def run_csis_pipeline(work_dir: Path) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, pd.DataFrame, float]:
    """Run osmget -> osmactive -> netfit (netclean) -> osmactive continuity -> netweld."""
    t0 = time.perf_counter()
    pbf = CASE_DIR / "source.osm.pbf"
    ways_pq = work_dir / "csis_ways.parquet"
    bar_raw_pq = work_dir / "csis_barriers_raw.parquet"
    ways_acc_pq = work_dir / "csis_ways_access.parquet"
    bar_pq = work_dir / "csis_barriers.parquet"
    ways_mot_pq = work_dir / "csis_ways_motor.parquet"
    links_nc_pq = work_dir / "csis_links_netclean.parquet"
    nodes_nc_pq = work_dir / "csis_nodes_motor.parquet"
    rep_json = work_dir / "csis_netclean_report.json"
    links_cont_pq = work_dir / "csis_links_cont.parquet"
    src_27700_pq = work_dir / "csis_ways_27700.parquet"
    tgt_27700_pq = work_dir / "csis_n3_27700_in.parquet"
    welded_pq = work_dir / "csis_n3_welded.parquet"

    subprocess.run([str(OSMGET_BIN), "extract", "-p", "all", "-i", str(pbf), "-o", str(ways_pq)], check=True)
    subprocess.run([str(OSMGET_BIN), "barriers", str(pbf), str(bar_raw_pq)], check=True)
    subprocess.run(
        [str(OSMACTIVE_BIN), "access", "--input", str(ways_pq), "--output", str(ways_acc_pq), "--rules", str(OSMACTIVE_RULES)],
        check=True,
    )
    subprocess.run(
        [str(OSMACTIVE_BIN), "barriers", "--input", str(bar_raw_pq), "--output", str(bar_pq), "--rules", str(OSMACTIVE_RULES)],
        check=True,
    )

    # Filter ways to motor-accessible before topological consolidation so footpaths don't split motor corridors
    w_acc = gpd.read_parquet(ways_acc_pq)
    w_mot = w_acc[w_acc["allow_car"].fillna(False).astype(bool)].copy()
    w_mot.to_parquet(ways_mot_pq, index=False)

    subprocess.run(
        [
            str(NETFIT_BIN), "network",
            "--centre", f"{CENTRE_LON},{CENTRE_LAT}",
            "--radius-m", str(RADIUS_M),
            "--osm", str(ways_mot_pq),
            "--out", str(links_nc_pq),
            "--barriers", str(bar_pq),
            "--nodes-out", str(nodes_nc_pq),
            "--access-rules", str(OSMACTIVE_RULES),
            "--report", str(rep_json),
        ],
        check=True,
    )
    subprocess.run(
        [str(OSMACTIVE_BIN), "continuity", "--input", str(links_nc_pq), "--output", str(links_cont_pq), "--rules", str(OSMACTIVE_RULES)],
        check=True,
    )

    # Prepare projected source & target for netweld transfer of directional lanes & posted maxspeed
    w_mot_27700 = w_mot.to_crs(27700)
    w_mot_27700["maxspeed_kmh"] = w_mot_27700["maxspeed"].apply(parse_maxspeed_kmh)
    w_mot_27700["dir_lanes"] = parse_dir_lanes(w_mot_27700)
    w_mot_27700[["osm_id", "maxspeed_kmh", "dir_lanes", "geometry"]].to_parquet(src_27700_pq, index=False)

    n3_in = gpd.read_parquet(links_cont_pq).to_crs(27700)
    n3_in.to_parquet(tgt_27700_pq, index=False)

    subprocess.run(
        [
            str(NETWELD_BIN), "join",
            str(src_27700_pq), str(tgt_27700_pq),
            "-o", str(welded_pq),
            "--distance", "12", "--angle", "25",
            "--transfer", "maxspeed_kmh:mean,dir_lanes:mean",
            "--prefix", "nw_",
        ],
        check=True,
    )
    elapsed_s = time.perf_counter() - t0

    n3 = gpd.read_parquet(welded_pq)
    if n3.geometry.name != "geometry":
        n3 = n3.rename_geometry("geometry")
    n3["segment_id"] = n3["link_id"].astype(str)
    n3["road_class"] = n3["highway"].map(HIGHWAY_TO_CLASS).fillna("minor")
    n3["name"] = n3["road_name"]
    n3["ref"] = n3["road_ref"]
    n3["dir_car"] = n3["dir_car"].map({0: "both", 1: "forward", -1: "backward"}).fillna("both")
    n3["allow_car"] = n3["allow_car"].fillna(True).astype(bool)
    n3["allow_heavy"] = n3["allow_heavy"].fillna(True).astype(bool)
    n3 = n3[n3["allow_car"]].copy().reset_index(drop=True)
    n3["length_m"] = n3.geometry.length

    nodes_bar = pd.read_parquet(nodes_nc_pq)
    closed_car = nodes_bar[~nodes_bar["allow_car"].fillna(True)]
    blocked_xy = None
    if len(closed_car) > 0:
        gbar = gpd.GeoSeries(gpd.points_from_xy(closed_car["lon"], closed_car["lat"]), crs=4326).to_crs(27700)
        blocked_xy = np.c_[gbar.x.values, gbar.y.values]

    # N3 keeps single weakly connected motor component; N4 severs car-blocking modal filters & prunes service/track ways
    u3, v3 = assign_endpoint_nodes(n3, blocked_xy=None)
    n3["u"], n3["v"] = u3, v3
    nd3 = np.unique(n3[["u", "v"]].values)
    idx3 = pd.Series(np.arange(len(nd3)), index=nd3)
    g3 = coo_matrix((np.ones(len(n3)), (idx3[n3["u"]].values, idx3[n3["v"]].values)), shape=(len(nd3), len(nd3)))
    _, lab3 = connected_components(g3, directed=True, connection="weak")
    lcc3 = set(nd3[lab3 == np.bincount(lab3).argmax()])
    n3 = n3[n3["u"].isin(lcc3) & n3["v"].isin(lcc3)].reset_index(drop=True)

    keep_n4 = ~n3["highway"].isin(["service", "track", "footway", "path", "cycleway", "pedestrian", "steps", "bridleway"])
    n4 = n3.loc[keep_n4].copy().reset_index(drop=True)
    u4, v4 = assign_endpoint_nodes(n4, blocked_xy=blocked_xy)
    n4["u"], n4["v"] = u4, v4
    # Retain largest weakly connected motor component on N4 after barrier severance
    nd4 = np.unique(n4[["u", "v"]].values)
    idx4 = pd.Series(np.arange(len(nd4)), index=nd4)
    g4 = coo_matrix((np.ones(len(n4)), (idx4[n4["u"]].values, idx4[n4["v"]].values)), shape=(len(nd4), len(nd4)))
    _, lab4 = connected_components(g4, directed=True, connection="weak")
    lcc4 = set(nd4[lab4 == np.bincount(lab4).argmax()])
    n4 = n4[n4["u"].isin(lcc4) & n4["v"].isin(lcc4)].reset_index(drop=True)

    rep_data = json.loads(rep_json.read_text())
    step_df = pd.DataFrame(rep_data["steps"])
    step_df["step_order"] = np.arange(1, len(step_df) + 1)
    step_df["input_ways"] = int(rep_data.get("ways", len(w_mot)))
    step_df["output_links"] = int(rep_data.get("links", len(n3)))
    step_df["blocked_barriers"] = int(len(closed_car))

    return n3, n4, step_df, elapsed_s


def network_stats(gdf: gpd.GeoDataFrame, variant_id: str, label: str, tools: str, lad_poly, build_s: float) -> list[dict]:
    """Compute topological and geometric metrics for both Full 15 km Cordon and Clipped Leeds LAD."""
    rows = []
    for scope_name, sub in [("full_15km_cordon", gdf), ("leeds_lad", gdf[gdf.intersects(lad_poly)].copy())]:
        if scope_name == "leeds_lad":
            sub["eff_len_m"] = sub.geometry.intersection(lad_poly).length
        else:
            sub["eff_len_m"] = sub["length_m"].values

        allow = sub[sub["allow_car"].fillna(False).astype(bool)]
        n_links = int(len(allow))
        net_km = float(allow["eff_len_m"].sum() / 1000.0)
        maj_km = float(allow.loc[allow["road_class"].isin(["motorway", "trunk", "primary"]), "eff_len_m"].sum() / 1000.0)
        mean_len_m = float(allow["length_m"].mean()) if n_links else 0.0

        dcar = allow["dir_car"].fillna("both").astype(str)
        n_dir_edges = int((dcar.isin(["both", "forward"])).sum() + (dcar.isin(["both", "backward"])).sum())

        nodes = np.unique(allow[["u", "v"]].values)
        n_nodes = int(len(nodes))
        if n_nodes > 0 and n_links > 0:
            idx = pd.Series(np.arange(n_nodes), index=nodes)
            g = coo_matrix((np.ones(n_links), (idx[allow["u"]].values, idx[allow["v"]].values)), shape=(n_nodes, n_nodes))
            n_comp, lab = connected_components(g, directed=True, connection="weak")
            top_lab = int(np.bincount(lab).argmax())
            lcc_node_set = set(nodes[lab == top_lab])
            in_lcc = allow["u"].isin(lcc_node_set) & allow["v"].isin(lcc_node_set)
            lcc_link_share = float(in_lcc.mean() * 100.0)
            lcc_km_share = float(allow.loc[in_lcc, "eff_len_m"].sum() / max(allow["eff_len_m"].sum(), 1e-9) * 100.0)
        else:
            n_comp, lcc_link_share, lcc_km_share = 0, 0.0, 0.0

        rows.append(
            {
                "variant_id": variant_id,
                "variant_label": label,
                "scope": scope_name,
                "tools": tools,
                "n_links": n_links,
                "n_directed_edges": n_dir_edges,
                "n_nodes": n_nodes,
                "network_km": net_km,
                "major_road_km": maj_km,
                "n_components": int(n_comp),
                "lcc_link_share_pct": lcc_link_share,
                "lcc_km_share_pct": lcc_km_share,
                "mean_link_length_m": mean_len_m,
                "build_time_s": round(build_s, 2),
            }
        )
    return rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    lad_poly = gpd.read_file(DATA_DIR / "leeds-lad.geojson").to_crs(27700).geometry.iloc[0]
    cx, cy = Transformer.from_crs(4326, 27700, always_xy=True).transform(CENTRE_LON, CENTRE_LAT)

    shas = {
        "aadtbench": git_sha(Path("/home/rbn/github/robinlovelace/aadtbench")),
        "osmget": git_sha(Path("/home/rbn/github/robinlovelace/osmget")),
        "osmactive": git_sha(Path("/home/rbn/github/atf-lab/osmactive")),
        "netfit": git_sha(Path("/home/rbn/github/atf-lab/netfit")),
        "netweld": git_sha(Path("/home/rbn/github/atfutures/netweld")),
    }
    print("Tool commits:", shas)

    # Load raw leeds-v1 network
    t0 = time.perf_counter()
    raw_all = gpd.read_parquet(CASE_DIR / "network.parquet")
    t_load = time.perf_counter() - t0

    # N0: Uncritical raw OSM (all 230,330 highway ways treated as motor-accessible)
    n0 = raw_all.copy()
    n0["allow_car"] = True
    n0["dir_car"] = n0["dir_car"].fillna("both")

    # N1: Raw OSM motor-accessible (`allow_car == True`)
    n1 = raw_all[raw_all["allow_car"].fillna(False).astype(bool)].copy().reset_index(drop=True)

    # N2: Minimal Python clean
    t2 = time.perf_counter()
    n2 = build_minimal_clean(n1, cx, cy)
    dt2 = time.perf_counter() - t2

    # N3 & N4: csis Rust toolchain
    n3, n4, step_df, dt_csis = run_csis_pipeline(WORK_DIR)

    # Save working full-domain assignment networks (gitignored)
    keep_common = [
        "segment_id", "u", "v", "length_m", "road_class", "highway", "name", "ref",
        "allow_car", "allow_heavy", "dir_car", "nw_maxspeed_kmh_mean", "nw_dir_lanes_mean", "geometry",
    ]
    for name, gdf in [
        ("N0_raw_unfiltered", n0),
        ("N1_raw_motor", n1),
        ("N2_minimal_python", n2),
        ("N3_csis_netclean", n3),
        ("N4_csis_strategic", n4),
    ]:
        cols = [c for c in keep_common if c in gdf.columns]
        gdf[cols].to_parquet(WORK_DIR / f"{name}.parquet", index=False)

    # Compute summary statistics
    summary_rows = []
    summary_rows.extend(
        network_stats(
            n0, "N0_raw_unfiltered", "Raw OSM (unfiltered all-ways)",
            f"pyosmium / aadtbench ({shas['aadtbench']}), no modal filter", lad_poly, t_load,
        )
    )
    summary_rows.extend(
        network_stats(
            n1, "N1_raw_motor", "Raw OSM (motor-accessible)",
            f"pyosmium / aadtbench ({shas['aadtbench']}), allow_car filter only", lad_poly, t_load,
        )
    )
    summary_rows.extend(
        network_stats(
            n2, "N2_minimal_python", "Minimal clean (Python LCC + trim + deg-2 weld)",
            "geopandas + scipy.sparse.csgraph + shapely.line_merge", lad_poly, dt2,
        )
    )
    csis_attr = (
        f"osmget ({shas['osmget']}) + osmactive ({shas['osmactive']}) + "
        f"netfit/netclean ({shas['netfit']}) + netweld ({shas['netweld']})"
    )
    summary_rows.extend(
        network_stats(n3, "N3_csis_netclean", "CSIS Rust clean (netclean + osmactive + netweld)", csis_attr, lad_poly, dt_csis)
    )
    summary_rows.extend(
        network_stats(
            n4, "N4_csis_strategic", "CSIS strategic clean (N3 + modal-filter severance + service prune)",
            csis_attr + " + barrier severance", lad_poly, dt_csis,
        )
    )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_parquet(OUT_DIR / "network_variants_summary.parquet", index=False)
    step_df.to_parquet(OUT_DIR / "netclean_step_log.parquet", index=False)

    print("\n=== Network Variants Summary ===")
    cols_show = [
        "variant_id", "scope", "n_links", "network_km", "major_road_km",
        "n_components", "lcc_link_share_pct", "mean_link_length_m", "build_time_s",
    ]
    print(summary_df[cols_show].to_string(index=False))


if __name__ == "__main__":
    main()

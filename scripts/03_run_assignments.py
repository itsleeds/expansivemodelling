#!/usr/bin/env python3
"""03_run_assignments.py — Run AequilibraE assignment matrix, ablation study, and validation.

Evaluates baseline (`aon`, `ue_bfw` from `aadtbench`) and improved (`ue_improved_phys`,
`ue_improved_cal`) AequilibraE traffic assignment configurations across the five network
cleaning variants (`N0_raw_unfiltered` .. `N4_csis_strategic`), validates against DfT AADF
counts on 3 km spatial cross-validation blocks, and writes compact caches to `data/processed/`.

Outputs (<10 MB total):
  - data/processed/validation_metrics.parquet
  - data/processed/category_breakdown.parquet
  - data/processed/site_predictions.parquet
  - data/processed/leeds_lad_aadt.parquet
"""
from __future__ import annotations

import contextlib
import os
import re
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Transformer
from scipy import stats
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from shapely import STRtree
from sklearn.linear_model import Ridge

os.environ.setdefault("AEQ_SHOW_PROGRESS", "FALSE")
from aequilibrae.matrix import AequilibraeMatrix
from aequilibrae.paths import Graph, TrafficAssignment, TrafficClass

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
OUT_DIR = DATA_DIR / "processed"
WORK_DIR = OUT_DIR / "_work"
CASE_DIR = Path("/home/rbn/github/robinlovelace/aadtbench/cases/leeds-v1")

CENTRE_LON, CENTRE_LAT = -1.5437941, 53.7974185
THREADS = 32

# --- Exact aadtbench baseline constants (from adapters/_base.py & adapters/aequilibrae.py) ---
CLASS_SPEED_KMH = {
    "motorway": 100.0, "trunk": 80.0, "primary": 50.0, "secondary": 45.0,
    "tertiary": 40.0, "minor": 30.0, "service": 15.0, "cycleway": 15.0,
    "footway": 5.0, "other": 30.0,
}
CLASS_CAPACITY_VPH = {
    "motorway": 2200.0, "trunk": 1800.0, "primary": 900.0, "secondary": 600.0,
    "tertiary": 400.0, "minor": 250.0, "service": 100.0, "cycleway": 300.0,
    "footway": 300.0, "other": 300.0,
}
DEFAULT_CAPACITY_VPH = 300.0
PEAK_HOUR_FACTOR = 0.10
MOTOR_CONDUCTANCE = {
    "motorway": 3.5, "trunk": 2.8, "primary": 1.8, "secondary": 1.4,
    "tertiary": 1.1, "minor": 0.6, "service": 0.2, "other": 0.6,
}

# --- Improved multi-class physical assignment constants ---
DEFAULT_DIR_LANES = {
    "motorway": 3.0, "motorway_link": 1.3,
    "trunk": 2.0, "trunk_link": 1.2,
    "primary": 1.25, "primary_link": 1.0,
    "secondary": 1.10, "secondary_link": 1.0,
    "tertiary": 1.0, "tertiary_link": 1.0,
    "unclassified": 0.9, "residential": 0.8, "living_street": 0.5,
    "connector": 0.9, "service": 0.4, "track": 0.3,
}
PER_LANE_CAP_VPH = {
    "motorway": 2200.0, "motorway_link": 1500.0,
    "trunk": 1800.0, "trunk_link": 1200.0,
    "primary": 1000.0, "primary_link": 850.0,
    "secondary": 900.0, "secondary_link": 750.0,
    "tertiary": 650.0, "tertiary_link": 550.0,
    "unclassified": 450.0, "residential": 280.0, "living_street": 150.0,
    "connector": 400.0, "service": 120.0, "track": 80.0,
}
DEFAULT_SPEED_KMH = {
    "motorway": 105.0, "motorway_link": 65.0,
    "trunk": 82.0, "trunk_link": 55.0,
    "primary": 50.0, "primary_link": 42.0,
    "secondary": 46.0, "secondary_link": 40.0,
    "tertiary": 39.0, "tertiary_link": 34.0,
    "unclassified": 32.0, "residential": 25.0, "living_street": 16.0,
    "connector": 28.0, "service": 14.0, "track": 10.0,
}
HIERARCHY_PENALTY = {
    "motorway": 1.00, "motorway_link": 1.04,
    "trunk": 1.00, "trunk_link": 1.05,
    "primary": 1.00, "primary_link": 1.05,
    "secondary": 1.02, "secondary_link": 1.06,
    "tertiary": 1.06, "tertiary_link": 1.10,
    "unclassified": 1.14, "residential": 1.32, "living_street": 1.65,
    "connector": 1.20, "service": 2.00, "track": 2.50,
}

# --- Crosswalk matching (adapted from aadtbench/benchmark/crosswalk.py) ---
CLASS_RANK = {
    "motorway": 9, "trunk": 8, "primary": 7, "secondary": 6, "tertiary": 5,
    "minor": 4, "service": 3, "cycleway": 2, "footway": 1, "other": 0,
}


def _norm_name(s: object) -> str:
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ""
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def _name_keys(row_name: object, row_ref: object) -> set[str]:
    keys = {_norm_name(row_name)}
    for part in re.split(r"[;,]", str(row_ref) if row_ref is not None else ""):
        keys.add(_norm_name(part))
    keys.discard("")
    return keys


def _travel_bearing(geom, pt, direction: str) -> float:
    d = shapely.line_locate_point(geom, pt)
    L = shapely.length(geom)
    a = shapely.line_interpolate_point(geom, max(d - 5.0, 0.0))
    b = shapely.line_interpolate_point(geom, min(d + 5.0, L))
    ang = float(np.degrees(np.arctan2(shapely.get_y(b) - shapely.get_y(a), shapely.get_x(b) - shapely.get_x(a))))
    return ang + 180.0 if direction == "backward" else ang


def build_crosswalk(network_gdf: gpd.GeoDataFrame, sites_df: pd.DataFrame, radius_m: float = 30.0) -> pd.DataFrame:
    """Match count points to motor-accessible network links with opposite-carriageway pairing."""
    net = network_gdf.reset_index(drop=True)
    geoms = net.geometry.values
    tree = STRtree(geoms)
    seg_ids = net["segment_id"].astype(str).values
    rank = net["road_class"].map(CLASS_RANK).fillna(0).astype(int).values
    n_col = net["name"] if "name" in net.columns else [None] * len(net)
    r_col = net["ref"] if "ref" in net.columns else [None] * len(net)
    name_keys = [_name_keys(n, r) for n, r in zip(n_col, r_col)]
    mode_ok = net["allow_car"].fillna(False).astype(bool).values
    dirs = net["dir_car"].fillna("both").astype(str).values

    rows = []
    for rec in sites_df.itertuples(index=False):
        pt = shapely.points(rec.x, rec.y)
        cand = tree.query(pt, predicate="dwithin", distance=radius_m)
        cand = cand[mode_ok[cand]]
        if len(cand) == 0:
            continue
        dist = shapely.distance(geoms[cand], pt)
        site_name = _norm_name(getattr(rec, "road_name", None))
        is_named = np.zeros(len(cand), dtype=bool)
        if site_name:
            is_named = np.array([site_name in name_keys[c] for c in cand], dtype=bool)
            if is_named.any():
                cand, dist, is_named = cand[is_named], dist[is_named], is_named[is_named]
        win = dist <= dist.min() + 10.0
        cand, dist, is_named = cand[win], dist[win], is_named[win]
        top = rank[cand] == rank[cand].max()
        cand, dist, is_named = cand[top], dist[top], is_named[top]
        tied = dist <= dist.min() + 1.0
        cand, dist, is_named = cand[tied], dist[tied], is_named[tied]
        lengths = shapely.length(geoms[cand])
        k = int(np.lexsort((dist, -lengths, -is_named.astype(int)))[0])
        i, d0 = int(cand[k]), float(dist[k])
        rows.append((str(rec.site_id), seg_ids[i], 1.0, round(d0, 2), "point"))

        # Check opposite dual carriageway partner within 50 m and bearing > 135 deg
        if dirs[i] != "both":
            pcand = tree.query(pt, predicate="dwithin", distance=50.0)
            pcand = pcand[(pcand != i) & mode_ok[pcand] & (dirs[pcand] != "both")]
            if len(pcand) > 0:
                if name_keys[i]:
                    pcand = pcand[np.array([bool(name_keys[i] & name_keys[c]) for c in pcand], dtype=bool)]
                else:
                    pcand = pcand[rank[pcand] == rank[i]]
            if len(pcand) > 0:
                b0 = _travel_bearing(geoms[i], pt, dirs[i])
                ang_ok = np.array(
                    [abs((_travel_bearing(geoms[c], pt, dirs[c]) - b0 + 180.0) % 360.0 - 180.0) > 135.0 for c in pcand],
                    dtype=bool,
                )
                pcand = pcand[ang_ok]
                if len(pcand) > 0:
                    j = int(pcand[np.argmin(shapely.distance(geoms[pcand], pt))])
                    rows.append(
                        (str(rec.site_id), seg_ids[j], 1.0, round(float(shapely.distance(geoms[j], pt)), 2), "pair")
                    )
    return pd.DataFrame(rows, columns=["site_id", "segment_id", "weight", "dist_m", "rule"])


# --- Baseline AequilibraE assignment helpers (exact aadtbench contract + nocliff option) ---
def build_baseline_edges_and_od(
    net: gpd.GeoDataFrame, zones: gpd.GeoDataFrame, od_motor: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    allow = net["allow_car"].fillna(False).astype(bool)
    seg = net.loc[allow, ["segment_id", "u", "v", "length_m", "road_class"]].copy()
    dcol = net.loc[allow, "dir_car"].fillna("both").astype(str)
    fwd = seg[dcol.isin(["both", "forward"])]
    bwd = seg[dcol.isin(["both", "backward"])].rename(columns={"u": "v", "v": "u"})
    edges = pd.concat([fwd, bwd[fwd.columns]], ignore_index=True)
    kmh = edges["road_class"].map(CLASS_SPEED_KMH).fillna(30.0)
    edges["time_s"] = edges["length_m"] / (kmh / 3.6)

    # LCC nodes & node coordinates
    ids = np.unique(edges[["u", "v"]].values)
    idx = pd.Series(np.arange(len(ids)), index=ids)
    g = coo_matrix((np.ones(len(edges)), (idx[edges["u"]].values, idx[edges["v"]].values)), shape=(len(ids), len(ids)))
    _, lab = connected_components(g, directed=True, connection="weak")
    lcc = ids[lab == np.bincount(lab).argmax()]

    first = shapely.get_point(net.geometry.values, 0)
    last = shapely.get_point(net.geometry.values, -1)
    node_xy = (
        pd.concat(
            [
                pd.DataFrame({"node": net["u"].values, "x": shapely.get_x(first), "y": shapely.get_y(first)}),
                pd.DataFrame({"node": net["v"].values, "x": shapely.get_x(last), "y": shapely.get_y(last)}),
            ]
        )
        .drop_duplicates("node")
        .set_index("node")
    )
    all_lcc_xy = node_xy.loc[lcc]

    z_vol = od_motor.groupby("o_zone")["trips"].sum().add(od_motor.groupby("d_zone")["trips"].sum(), fill_value=0.0)
    z_sorted = zones.assign(_vol=zones["zone_id"].map(z_vol).fillna(0.0)).sort_values("_vol", ascending=False)
    top_z, rest_z = z_sorted.iloc[:300], z_sorted.iloc[300:]

    def _snap(z_sub: gpd.GeoDataFrame, ndf: pd.DataFrame, gn: int) -> pd.DataFrame:
        if z_sub.empty:
            return pd.DataFrame(columns=["zone_id", "node", "w"])
        tree = cKDTree(ndf[["x", "y"]].values)
        nids = ndf.index.values
        b = z_sub.bounds
        zids = z_sub["zone_id"].values
        offs = np.linspace(1.0 / (2 * gn), 1.0 - 1.0 / (2 * gn), gn)
        ox, oy = np.meshgrid(offs, offs)
        ox, oy = ox.ravel(), oy.ravel()
        xs = (b["minx"].values[:, None] + (b["maxx"].values - b["minx"].values)[:, None] * ox[None, :]).ravel()
        ys = (b["miny"].values[:, None] + (b["maxy"].values - b["miny"].values)[:, None] * oy[None, :]).ravel()
        d, ii = tree.query(np.c_[xs, ys])
        ok = d <= 1000.0
        z_rep = np.repeat(zids, len(ox))
        zc = pd.DataFrame({"zone_id": z_rep[ok], "node": nids[ii[ok]], "w": 1.0}).groupby(["zone_id", "node"], as_index=False)["w"].sum()
        zc["w"] /= zc.groupby("zone_id")["w"].transform("sum")
        return zc

    coll = edges[edges["road_class"].isin({"motorway", "trunk", "primary", "secondary", "tertiary"})]
    pref_ids = np.intersect1d(lcc, np.unique(coll[["u", "v"]].values))
    if len(pref_ids) >= 20:
        zc = pd.concat([_snap(top_z, node_xy.loc[pref_ids], 2), _snap(rest_z, node_xy.loc[pref_ids], 1)], ignore_index=True)
        missed = z_sorted[~z_sorted["zone_id"].isin(set(zc["zone_id"]))]
        if not missed.empty:
            zc = pd.concat([zc, _snap(missed, all_lcc_xy, 1)], ignore_index=True)
    else:
        zc = pd.concat([_snap(top_z, all_lcc_xy, 2), _snap(rest_z, all_lcc_xy, 1)], ignore_index=True)

    m = od_motor.merge(zc.rename(columns={"zone_id": "o_zone", "node": "o_node", "w": "ow"}), on="o_zone")
    m = m.merge(zc.rename(columns={"zone_id": "d_zone", "node": "d_node", "w": "dw"}), on="d_zone")
    m["trips"] = m["trips"] * m["ow"] * m["dw"]
    m = m[m["o_node"] != m["d_node"]]
    od_nodes = m.groupby(["o_node", "d_node"], as_index=False)["trips"].sum()
    return edges.reset_index(drop=True), od_nodes


def run_single_class_aeq(
    edges: pd.DataFrame,
    od_nodes: pd.DataFrame,
    algorithm: str = "all-or-nothing",
    max_iter: int = 20,
    rgap: float = 1e-3,
    use_cap_cliff: bool = True,
    custom_cap: np.ndarray | None = None,
    custom_time: np.ndarray | None = None,
) -> tuple[pd.DataFrame, float, float, float]:
    t0 = time.perf_counter()
    ue = algorithm != "all-or-nothing"
    scale = PEAK_HOUR_FACTOR if ue else 1.0
    raw = np.unique(edges[["u", "v"]].values.astype(np.int64))
    node_map = pd.Series(np.arange(1, len(raw) + 1, dtype=np.int64), index=raw)
    cap = custom_cap if custom_cap is not None else edges["road_class"].map(CLASS_CAPACITY_VPH).fillna(DEFAULT_CAPACITY_VPH).values
    ff_t = custom_time if custom_time is not None else edges["time_s"].values
    cap_km = float((cap * edges["length_m"].values).sum()) / 1000.0
    trip_sum = float(od_nodes["trips"].values.sum())
    if ue and use_cap_cliff and (trip_sum * scale > 0.09 * cap_km) and trip_sum > 0:
        scale = 0.005 * cap_km / trip_sum

    net_df = pd.DataFrame(
        {
            "link_id": np.arange(1, len(edges) + 1, dtype=np.int64),
            "a_node": node_map.loc[edges["u"].values].values,
            "b_node": node_map.loc[edges["v"].values].values,
            "direction": np.int8(1),
            "distance": edges["length_m"].values.astype(float),
            "free_flow_time": np.maximum(ff_t.astype(float), 1e-3),
            "capacity": np.maximum(cap.astype(float), 50.0),
        }
    )
    o = node_map.loc[od_nodes["o_node"].values].values
    d = node_map.loc[od_nodes["d_node"].values].values
    used = np.unique(np.concatenate([o, d]))
    pos = pd.Series(np.arange(len(used)), index=used)
    arr = np.zeros((len(used), len(used)))
    np.add.at(arr, (pos.loc[o].values, pos.loc[d].values), od_nodes["trips"].values * scale)

    mat = AequilibraeMatrix()
    mat.create_empty(zones=len(used), matrix_names=["demand"], index_names=["node_id"], memory_only=True)
    mat.index[:] = used
    mat.matrix["demand"][:, :] = arr
    mat.computational_view(["demand"])

    graph = Graph()
    graph.network = net_df
    graph.mode = "c"
    graph.prepare_graph(centroids=used)
    graph.set_blocked_centroid_flows(False)
    graph.set_graph("free_flow_time")

    tc = TrafficClass("cars", graph, mat)
    tc.set_pce(1.0)
    assig = TrafficAssignment()
    assig.set_classes([tc])
    assig.set_vdf("BPR")
    assig.set_vdf_parameters({"alpha": 0.15, "beta": 4.0})
    assig.set_capacity_field("capacity")
    assig.set_time_field("free_flow_time")
    assig.max_iter = max_iter
    assig.rgap_target = rgap
    assig.set_algorithm(algorithm)
    assig.set_cores(THREADS)
    with open(os.devnull, "w") as sink, contextlib.redirect_stderr(sink):
        assig.execute(log_specification=False)

    rep = assig.report()
    final_rgap = float(rep["rgap"].iloc[-1]) if (ue and len(rep) > 0 and "rgap" in rep.columns) else 0.0
    edge_flow = np.maximum(assig.results()["demand_tot"].reindex(np.arange(1, len(edges) + 1)).fillna(0.0).values, 0.0) / scale
    seg_flow = pd.DataFrame({"segment_id": edges["segment_id"].astype(str).values, "flow": edge_flow}).groupby("segment_id", as_index=False)["flow"].sum()
    dt = time.perf_counter() - t0
    return seg_flow, dt, scale, final_rgap


def modulate_aadtbench(net: gpd.GeoDataFrame, seg_flow: pd.DataFrame) -> pd.DataFrame:
    """Exact aadtbench motor modulation: ((flow + 0.03 * median_pos) * cond^0.7)^0.5."""
    rc_map = net.set_index("segment_id")["road_class"]
    out = seg_flow.copy()
    fl = np.maximum(out["flow"].to_numpy(float), 0.0)
    if (fl > 0).any():
        med = float(np.median(fl[fl > 0]))
        cond = out["segment_id"].map(rc_map).map(MOTOR_CONDUCTANCE).fillna(1.0).to_numpy(float) ** 0.7
        out["flow"] = ((fl + 0.03 * med) * cond) ** 0.50
    return out


# --- Improved multi-class SRN-cordon BFW UE ---
def enrich_network_attributes(net: gpd.GeoDataFrame, pois: pd.DataFrame, pop: pd.DataFrame) -> gpd.GeoDataFrame:
    out = net.copy()
    hw = (out["highway"] if "highway" in out.columns else out["road_class"]).fillna("unclassified").astype(str)
    def_lanes = hw.map(DEFAULT_DIR_LANES).fillna(1.0)
    dir_lanes = out["nw_dir_lanes_mean"].fillna(def_lanes) if "nw_dir_lanes_mean" in out.columns else def_lanes
    out["dir_lanes"] = dir_lanes.astype(float)
    out["cap_vph"] = (out["dir_lanes"] * hw.map(PER_LANE_CAP_VPH).fillna(450.0)).astype(float)

    def_spd = hw.map(DEFAULT_SPEED_KMH).fillna(30.0).to_numpy(float)
    if "nw_maxspeed_kmh_mean" in out.columns:
        spd_raw = out["nw_maxspeed_kmh_mean"].to_numpy(float)
        is_mw = hw.isin(["motorway", "trunk"]).to_numpy(bool)
        spd = np.where(np.isfinite(spd_raw), spd_raw * np.where(is_mw, 0.96, 0.85), def_spd)
    else:
        spd = def_spd
    out["speed_kmh"] = np.clip(spd, 10.0, 115.0)
    out["hier_pen"] = hw.map(HIERARCHY_PENALTY).fillna(1.14).astype(float)

    # Spatial kernel density of WorldPop + OSM POIs (open T1 covariates in leeds-v1/inputs/)
    mid = shapely.line_interpolate_point(out.geometry.values, 0.5, normalized=True)
    mx, my = shapely.get_x(mid), shapely.get_y(mid)
    act_xy = np.vstack([pois[["x", "y"]].values, pop[["x", "y"]].values])
    act_w = np.r_[pois["weight"].values, pop["weight"].values * 0.25]
    cell_m = 250.0
    df_act = (
        pd.DataFrame({"gx": np.floor(act_xy[:, 0] / cell_m).astype(int), "gy": np.floor(act_xy[:, 1] / cell_m).astype(int), "w": act_w})
        .groupby(["gx", "gy"], as_index=False)["w"]
        .sum()
    )
    l_cells, inv_c = np.unique(np.c_[np.floor(mx / cell_m).astype(int), np.floor(my / cell_m).astype(int)], axis=0, return_inverse=True)
    pairs = cKDTree((l_cells + 0.5) * cell_m).sparse_distance_matrix(
        cKDTree((df_act[["gx", "gy"]].values + 0.5) * cell_m), 1500.0, output_type="ndarray"
    )
    attr = np.bincount(pairs["i"], df_act["w"].values[pairs["j"]] * np.exp(-pairs["v"] / 500.0), minlength=len(l_cells))[inv_c.ravel()]
    med_a = float(np.median(attr[attr > 0])) if (attr > 0).any() else 1.0
    out["log_attr"] = np.log1p(attr / max(med_a, 1.0))
    return out


def run_improved_multiclass_ue(
    net_enr: gpd.GeoDataFrame,
    zones: gpd.GeoDataFrame,
    od_raw: pd.DataFrame,
    cx: float,
    cy: float,
) -> tuple[pd.DataFrame, float, float, int]:
    """Run multi-class (`car` PCE=1, `heavy` PCE=2) BFW UE with SRN cordon gateways and `netweld` capacities/speeds."""
    t0 = time.perf_counter()
    allow = net_enr["allow_car"].fillna(False).astype(bool)
    sub = net_enr.loc[allow].copy().reset_index(drop=True)
    hw_col = (sub["highway"] if "highway" in sub.columns else sub["road_class"]).fillna("unclassified").astype(str)
    sub["_hw"] = hw_col

    dcar = sub["dir_car"].fillna("both").astype(str)
    fwd = sub[dcar.isin(["both", "forward"])].copy()
    bwd = sub[dcar.isin(["both", "backward"])].copy().rename(columns={"u": "v", "v": "u"})
    edges = pd.concat([fwd, bwd[fwd.columns]], ignore_index=True)

    # Pure physical free-flow travel time from netweld posted maxspeed / road hierarchy speed
    # (HGVs capped at UK 90 km/h speed limit, slight residential turn impedance for heavy goods only)
    edges["ff_time_s"] = edges["length_m"].values / (edges["speed_kmh"].values / 3.6)
    edges["ff_time_heavy_s"] = (
        edges["length_m"].values / (np.minimum(edges["speed_kmh"].values, 90.0) / 3.6)
    ) * np.where(edges["_hw"].isin(["residential", "living_street"]), 1.6, 1.0)

    # LCC filter
    e_nodes = np.unique(edges[["u", "v"]].values)
    idx_map = pd.Series(np.arange(len(e_nodes)), index=e_nodes)
    g = coo_matrix((np.ones(len(edges)), (idx_map.loc[edges["u"]].values, idx_map.loc[edges["v"]].values)), shape=(len(e_nodes), len(e_nodes)))
    _, labels = connected_components(g, directed=True, connection="weak")
    lcc_nodes = set(e_nodes[labels == np.bincount(labels).argmax()])
    edges = edges[edges["u"].isin(lcc_nodes) & edges["v"].isin(lcc_nodes)].reset_index(drop=True)
    lcc_arr = np.array(sorted(lcc_nodes))

    pts_u = shapely.get_point(sub.geometry.values, 0)
    pts_v = shapely.get_point(sub.geometry.values, -1)
    node_xy = (
        pd.concat(
            [
                pd.DataFrame({"node": sub["u"].values, "x": shapely.get_x(pts_u), "y": shapely.get_y(pts_u)}),
                pd.DataFrame({"node": sub["v"].values, "x": shapely.get_x(pts_v), "y": shapely.get_y(pts_v)}),
            ]
        )
        .drop_duplicates("node")
        .set_index("node")
        .loc[lcc_arr]
    )

    # Collector-first zone connectors (matching baseline spatial resolution: top 300 zones 2x2, rest 1x1)
    od_car = od_raw[od_raw["mode"] == "car"][["o_zone", "d_zone", "trips"]]
    od_hvy = od_raw[od_raw["mode"] == "heavy"][["o_zone", "d_zone", "trips"]]
    z_vol = od_car.groupby("o_zone")["trips"].sum().add(od_car.groupby("d_zone")["trips"].sum(), fill_value=0.0)
    z_sorted = zones.assign(_vol=zones["zone_id"].map(z_vol).fillna(0.0)).sort_values("_vol", ascending=False)
    top_z, rest_z = z_sorted.iloc[:300], z_sorted.iloc[300:]

    def _snap_imp(z_sub: gpd.GeoDataFrame, ndf: pd.DataFrame, gn: int) -> pd.DataFrame:
        if z_sub.empty:
            return pd.DataFrame(columns=["zone_id", "node", "w"])
        tree = cKDTree(ndf[["x", "y"]].values)
        nids = ndf.index.values
        b = z_sub.bounds
        zids = z_sub["zone_id"].values
        offs = np.linspace(1.0 / (2 * gn), 1.0 - 1.0 / (2 * gn), gn)
        ox, oy = np.meshgrid(offs, offs)
        ox, oy = ox.ravel(), oy.ravel()
        xs = (b["minx"].values[:, None] + (b["maxx"].values - b["minx"].values)[:, None] * ox[None, :]).ravel()
        ys = (b["miny"].values[:, None] + (b["maxy"].values - b["miny"].values)[:, None] * oy[None, :]).ravel()
        d, ii = tree.query(np.c_[xs, ys])
        ok = d <= 1000.0
        z_rep = np.repeat(zids, len(ox))
        zc_sub = pd.DataFrame({"zone_id": z_rep[ok], "node": nids[ii[ok]], "w": 1.0}).groupby(["zone_id", "node"], as_index=False)["w"].sum()
        zc_sub["w"] /= zc_sub.groupby("zone_id")["w"].transform("sum")
        return zc_sub

    coll_e = edges[edges["road_class"].isin({"primary", "secondary", "tertiary"}) | edges["_hw"].isin({"unclassified"})]
    pref_ids = np.intersect1d(lcc_arr, np.unique(coll_e[["u", "v"]].values))
    if len(pref_ids) >= 20:
        zc = pd.concat([_snap_imp(top_z, node_xy.loc[pref_ids], 2), _snap_imp(rest_z, node_xy.loc[pref_ids], 1)], ignore_index=True)
        missed = z_sorted[~z_sorted["zone_id"].isin(set(zc["zone_id"]))]
        if not missed.empty:
            zc = pd.concat([zc, _snap_imp(missed, node_xy, 1)], ignore_index=True)
    else:
        zc = pd.concat([_snap_imp(top_z, node_xy, 2), _snap_imp(rest_z, node_xy, 1)], ignore_index=True)

    def _map_od(df_z: pd.DataFrame) -> pd.DataFrame:
        m = df_z.merge(zc.rename(columns={"zone_id": "o_zone", "node": "o_node", "w": "ow"}), on="o_zone")
        m = m.merge(zc.rename(columns={"zone_id": "d_zone", "node": "d_node", "w": "dw"}), on="d_zone")
        m["trips"] = m["trips"] * m["ow"] * m["dw"]
        m = m[m["o_node"] != m["d_node"]]
        return m.groupby(["o_node", "d_node"], as_index=False)["trips"].sum()

    nod_int_car = _map_od(od_car)
    nod_int_hvy = _map_od(od_hvy)

    # SRN perimeter gateways (`motorway` and `trunk` portals at r >= 14.2 km, clustered in 18 deg sectors)
    node_xy["r"] = np.hypot(node_xy["x"].values - cx, node_xy["y"].values - cy)
    srn_e = edges[edges["road_class"].isin({"motorway", "trunk"})].copy()
    srn_e["ru"] = node_xy.loc[srn_e["u"].values, "r"].values
    srn_e["rv"] = node_xy.loc[srn_e["v"].values, "r"].values
    outer = pd.concat(
        [
            srn_e[srn_e["ru"] >= 14200.0][["u", "cap_vph", "road_class", "ru"]].rename(columns={"u": "node", "ru": "r"}),
            srn_e[srn_e["rv"] >= 14200.0][["v", "cap_vph", "road_class", "rv"]].rename(columns={"v": "node", "rv": "r"}),
        ],
        ignore_index=True,
    )
    outer["x"] = node_xy.loc[outer["node"].values, "x"].values
    outer["y"] = node_xy.loc[outer["node"].values, "y"].values
    outer["theta_bin"] = np.floor((np.degrees(np.arctan2(outer["y"] - cy, outer["x"] - cx)) + 180.0) / 18.0).astype(int)
    gw = outer.sort_values(["cap_vph", "r"], ascending=[False, False]).groupby("theta_bin", as_index=False).first()

    ang = np.arctan2(gw["y"].values - cy, gw["x"].values - cx)
    gx_arr, gy_arr, gcap, gnodes = gw["x"].values, gw["y"].values, gw["cap_vph"].values, gw["node"].values
    d_gg = np.hypot(gx_arr[:, None] - gx_arr[None, :], gy_arr[:, None] - gy_arr[None, :])
    dang = np.abs((ang[:, None] - ang[None, :] + np.pi) % (2 * np.pi) - np.pi)
    is_mw = (gw["road_class"] == "motorway").values
    mw_w = np.where(is_mw, 1.0, 0.20)
    valid_gg = (d_gg >= 14000.0) & (dang >= np.radians(80.0))
    w_gg = np.where(valid_gg, (gcap * mw_w)[:, None] * (gcap * mw_w)[None, :] * np.exp(-d_gg / 40000.0), 0.0)
    r_sum = np.maximum(w_gg.sum(axis=1, keepdims=True), 1.0)
    thru_mat = (gcap[:, None] * np.where(is_mw[:, None], 5.2, 0.85)) * (w_gg / r_sum)
    gi, gj = np.nonzero(thru_mat >= 1.0)
    nod_thru = pd.DataFrame({"o_node": gnodes[gi], "d_node": gnodes[gj], "trips": thru_mat[gi, gj]})

    z_cent = zones.assign(zx=zones.geometry.centroid.x, zy=zones.geometry.centroid.y)
    z_attr = z_cent["attractors"].values + 0.15 * z_cent["population"].values
    d_gz = np.hypot(gx_arr[:, None] - z_cent["zx"].values[None, :], gy_arr[:, None] - z_cent["zy"].values[None, :])
    w_gz = np.where(d_gz >= 4000.0, z_attr[None, :] * np.exp(-d_gz / 14000.0), 0.0)
    rz_sum = np.maximum(w_gz.sum(axis=1, keepdims=True), 1.0)
    ext_one_way = (gcap[:, None] * np.where(is_mw[:, None], 2.4, 0.70)) * (w_gz / rz_sum)
    gi2, zj2 = np.nonzero(ext_one_way >= 0.5)
    ext_in = pd.DataFrame({"g_node": gnodes[gi2], "zone_id": z_cent["zone_id"].values[zj2], "trips": ext_one_way[gi2, zj2]}).merge(zc, on="zone_id")
    ext_in["trips"] *= ext_in["w"]
    ext_in = ext_in[ext_in["g_node"] != ext_in["node"]]
    nod_ext = pd.concat(
        [
            ext_in[["g_node", "node", "trips"]].rename(columns={"g_node": "o_node", "node": "d_node"}),
            ext_in[["node", "g_node", "trips"]].rename(columns={"node": "o_node", "g_node": "d_node"}),
        ],
        ignore_index=True,
    )

    # Physical daily demand (~1.5M daily motor trips across the 15 km cordon, ~9.5% heavy vehicle share)
    car_all = (
        pd.concat(
            [
                nod_int_car.assign(trips=nod_int_car["trips"] * 0.62),
                nod_ext.assign(trips=nod_ext["trips"] * 0.78),
                nod_thru.assign(trips=nod_thru["trips"] * 0.78),
            ],
            ignore_index=True,
        )
        .groupby(["o_node", "d_node"], as_index=False)["trips"]
        .sum()
    )
    hvy_all = (
        pd.concat(
            [
                nod_int_hvy.assign(trips=nod_int_hvy["trips"] * 0.058),
                nod_ext.assign(trips=nod_ext["trips"] * 0.13),
                nod_thru.assign(trips=nod_thru["trips"] * 0.19),
            ],
            ignore_index=True,
        )
        .groupby(["o_node", "d_node"], as_index=False)["trips"]
        .sum()
    )

    phf = 0.085
    raw_nodes = np.unique(edges[["u", "v"]].values.astype(np.int64))
    node_map = pd.Series(np.arange(1, len(raw_nodes) + 1, dtype=np.int64), index=raw_nodes)
    used_orig = np.unique(np.r_[car_all["o_node"].values, car_all["d_node"].values, hvy_all["o_node"].values, hvy_all["d_node"].values])
    used = np.unique(node_map.loc[used_orig].values)
    pos = pd.Series(np.arange(len(used)), index=used)

    def _build_mat(df_od: pd.DataFrame, cname: str) -> AequilibraeMatrix:
        m = AequilibraeMatrix()
        m.create_empty(zones=len(used), matrix_names=[cname], index_names=["node_id"], memory_only=True)
        m.index[:] = used
        arr = np.zeros((len(used), len(used)))
        np.add.at(
            arr,
            (pos.loc[node_map.loc[df_od["o_node"].values].values].values, pos.loc[node_map.loc[df_od["d_node"].values].values].values),
            df_od["trips"].values * phf,
        )
        m.matrix[cname][:, :] = arr
        m.computational_view([cname])
        return m

    mat_car = _build_mat(car_all, "car")
    mat_hvy = _build_mat(hvy_all, "heavy")

    net_aeq = pd.DataFrame(
        {
            "link_id": np.arange(1, len(edges) + 1, dtype=np.int64),
            "a_node": node_map.loc[edges["u"].values].values,
            "b_node": node_map.loc[edges["v"].values].values,
            "direction": np.int8(1),
            "distance": edges["length_m"].values.astype(float),
            "ff_time_car": np.maximum(edges["ff_time_s"].values.astype(float), 1e-3),
            "ff_time_heavy": np.maximum(edges["ff_time_heavy_s"].values.astype(float), 1e-3),
            "capacity": np.maximum(edges["cap_vph"].values.astype(float), 50.0),
        }
    )

    g_car = Graph()
    g_car.network = net_aeq[["link_id", "a_node", "b_node", "direction", "distance", "ff_time_car", "capacity"]].rename(
        columns={"ff_time_car": "free_flow_time"}
    )
    g_car.mode = "c"
    g_car.prepare_graph(centroids=used)
    g_car.set_blocked_centroid_flows(False)
    g_car.set_graph("free_flow_time")

    hvy_ok = edges["allow_heavy"].fillna(True).astype(bool).values if "allow_heavy" in edges.columns else np.ones(len(edges), dtype=bool)
    g_hvy = Graph()
    g_hvy.network = net_aeq[["link_id", "a_node", "b_node", "direction", "distance", "capacity"]].assign(
        free_flow_time=np.where(hvy_ok, net_aeq["ff_time_heavy"].values, net_aeq["ff_time_heavy"].values * 50.0)
    )
    g_hvy.mode = "h"
    g_hvy.prepare_graph(centroids=used)
    g_hvy.set_blocked_centroid_flows(False)
    g_hvy.set_graph("free_flow_time")

    tc_car = TrafficClass("car", g_car, mat_car)
    tc_car.set_pce(1.0)
    tc_hvy = TrafficClass("heavy", g_hvy, mat_hvy)
    tc_hvy.set_pce(2.0)

    assig = TrafficAssignment()
    assig.set_classes([tc_car, tc_hvy])
    assig.set_vdf("BPR")
    assig.set_vdf_parameters({"alpha": 0.15, "beta": 4.0})
    assig.set_capacity_field("capacity")
    assig.set_time_field("free_flow_time")
    assig.max_iter = 25
    assig.rgap_target = 5e-4
    assig.set_algorithm("bfw")
    assig.set_cores(THREADS)
    with open(os.devnull, "w") as sink, contextlib.redirect_stderr(sink):
        assig.execute(log_specification=False)

    rep = assig.report()
    final_rgap = float(rep["rgap"].iloc[-1]) if (len(rep) > 0 and "rgap" in rep.columns) else 0.0
    res = assig.results().reindex(np.arange(1, len(edges) + 1)).fillna(0.0)
    pce_hr = np.maximum(res["car_tot"].values, 0.0) + 2.0 * np.maximum(res["heavy_tot"].values, 0.0)
    vc_edge = pce_hr / np.maximum(edges["cap_vph"].values, 50.0)
    flow_car = np.maximum(res["car_tot"].values, 0.0) / phf
    flow_hvy = np.maximum(res["heavy_tot"].values, 0.0) / phf

    seg_pred = (
        pd.DataFrame(
            {
                "segment_id": edges["segment_id"].astype(str).values,
                "flow": flow_car + flow_hvy,
                "car": flow_car,
                "heavy": flow_hvy,
                "vc_ratio": vc_edge,
            }
        )
        .groupby("segment_id", as_index=False)
        .agg({"flow": "sum", "car": "sum", "heavy": "sum", "vc_ratio": "max"})
    )
    dt = time.perf_counter() - t0
    return seg_pred, dt, final_rgap, int(len(gw))


# --- Spatially blocked calibration & scoring ---
RC_CATS = ["motorway", "trunk", "primary", "secondary", "tertiary", "minor"]


def _design_matrix(df_in: pd.DataFrame, flow_col: str = "flow") -> np.ndarray:
    fl = np.maximum(df_in[flow_col].to_numpy(float), 0.0)
    X = pd.DataFrame(
        {
            "log_flow": np.log(np.maximum(fl, 50.0)),
            "zero_flow": (fl < 50.0).astype(float),
            "log_cap2w": np.log(np.maximum(df_in["cap_2w"].to_numpy(float), 100.0)),
        }
    )
    rc = df_in["road_class"].astype(str).to_numpy()
    for c in RC_CATS[:-1]:
        X[f"rc_{c}"] = (rc == c).astype(float)
    return X.to_numpy(float)


def fit_oof_and_network_calibration(
    site_tbl: pd.DataFrame,
    net_enr: gpd.GeoDataFrame,
    seg_pred: pd.DataFrame,
    target_col: str = "obs_motor",
    flow_col: str = "flow",
) -> tuple[np.ndarray, np.ndarray]:
    """5-fold spatial OOF Ridge log-linear calibration with Duan smearing + full-fit link predictions."""
    oof_pred = np.zeros(len(site_tbl), dtype=float)
    y_all = np.log(np.maximum(site_tbl[target_col].to_numpy(float), 1.0))
    X_all = _design_matrix(site_tbl, flow_col=flow_col)

    for f in sorted(site_tbl["fold"].unique()):
        te = (site_tbl["fold"] == f).to_numpy(bool)
        tr = ~te
        reg = Ridge(alpha=1.0).fit(X_all[tr], y_all[tr])
        smear = float(np.mean(np.exp(y_all[tr] - reg.predict(X_all[tr]))))
        oof_pred[te] = np.exp(reg.predict(X_all[te])) * smear

    reg_full = Ridge(alpha=1.0).fit(X_all, y_all)
    smear_full = float(np.mean(np.exp(y_all - reg_full.predict(X_all))))

    net_fl = (
        net_enr[["segment_id"]]
        .assign(segment_id=net_enr["segment_id"].astype(str))
        .merge(seg_pred[["segment_id", flow_col]], on="segment_id", how="left")[flow_col]
        .fillna(0.0)
        .to_numpy(float)
    )
    is_1way = net_enr["dir_car"].fillna("both").astype(str).isin(["forward", "backward"]).to_numpy(bool)
    net_df_2w = pd.DataFrame(
        {
            flow_col: np.where(is_1way, net_fl * 2.0, net_fl),
            "cap_2w": net_enr["cap_vph"].to_numpy(float) * 2.0,
            "speed_kmh": net_enr["speed_kmh"].to_numpy(float),
            "log_attr": net_enr["log_attr"].to_numpy(float),
            "road_class": net_enr["road_class"].astype(str).to_numpy(),
        }
    )
    net_cal = np.exp(reg_full.predict(_design_matrix(net_df_2w, flow_col=flow_col))) * smear_full * np.where(is_1way, 0.5, 1.0)
    net_cal = np.where(net_enr["allow_car"].fillna(False).to_numpy(bool), net_cal, 0.0)
    return oof_pred, net_cal


def evaluate_predictions(
    sites_eval: pd.DataFrame,
    cw: pd.DataFrame,
    net_enr: gpd.GeoDataFrame,
    seg_flow: pd.DataFrame,
    lad_in_mask: np.ndarray,
    lad_len_m: np.ndarray,
    scaling_mode: str = "oof_vol",
) -> tuple[dict, pd.DataFrame, np.ndarray]:
    """Compute full suite of site validation metrics + Leeds LAD MVKM."""
    s = sites_eval.copy()
    s["site_id"] = s["site_id"].astype(str)
    pr = seg_flow[["segment_id", "flow"]].copy()
    pr["segment_id"] = pr["segment_id"].astype(str)

    xw = cw.merge(pr, on="segment_id", how="left")
    xw["wflow"] = xw["flow"].fillna(0.0).clip(lower=0.0) * xw["weight"].astype(float)
    agg = xw.groupby("site_id", as_index=False)["wflow"].sum().rename(columns={"wflow": "flow"})

    prim = cw[cw["rule"] != "pair"].merge(
        net_enr[["segment_id", "road_class", "cap_vph", "speed_kmh", "log_attr"]].assign(segment_id=net_enr["segment_id"].astype(str)),
        on="segment_id", how="left",
    )
    cap_2w = (
        cw.merge(net_enr[["segment_id", "cap_vph"]].assign(segment_id=net_enr["segment_id"].astype(str)), on="segment_id", how="left")
        .groupby("site_id")["cap_vph"]
        .sum()
    )
    tbl = s.merge(agg, on="site_id", how="inner")
    tbl["road_class"] = tbl["site_id"].map(prim.set_index("site_id")["road_class"]).fillna("minor")
    tbl["cap_2w"] = tbl["site_id"].map(cap_2w).fillna(900.0)
    tbl["speed_kmh"] = tbl["site_id"].map(prim.set_index("site_id")["speed_kmh"]).fillna(30.0)
    tbl["log_attr"] = tbl["site_id"].map(prim.set_index("site_id")["log_attr"]).fillna(0.0)

    net_raw_f = (
        net_enr[["segment_id"]]
        .assign(segment_id=net_enr["segment_id"].astype(str))
        .merge(pr, on="segment_id", how="left")["flow"]
        .fillna(0.0)
        .to_numpy(float)
    )

    obs = tbl["obs_motor"].to_numpy(float)
    raw_f = tbl["flow"].to_numpy(float)

    if scaling_mode == "none":
        pred = raw_f.copy()
        net_f_scaled = net_raw_f.copy()
        k_mean = 1.0
    elif scaling_mode == "oof_vol":
        pred = np.zeros(len(tbl), dtype=float)
        ks = []
        for f in sorted(tbl["fold"].unique()):
            te = (tbl["fold"] == f).to_numpy(bool)
            k = float(obs[~te].sum() / max(raw_f[~te].sum(), 1.0))
            pred[te] = k * raw_f[te]
            ks.append(k)
        k_mean = float(np.mean(ks))
        net_f_scaled = net_raw_f * k_mean
    elif scaling_mode == "oof_geom":
        pred = np.zeros(len(tbl), dtype=float)
        ks = []
        for f in sorted(tbl["fold"].unique()):
            te = (tbl["fold"] == f).to_numpy(bool)
            ok_tr = (~te) & (raw_f > 0) & (obs > 0)
            k = float(np.exp(np.mean(np.log(obs[ok_tr] / raw_f[ok_tr])))) if ok_tr.any() else 1.0
            pred[te] = k * raw_f[te]
            ks.append(k)
        k_mean = float(np.mean(ks))
        net_f_scaled = net_raw_f * k_mean
    elif scaling_mode == "oof_csis_cal":
        pred, net_f_scaled = fit_oof_and_network_calibration(tbl, net_enr, pr, target_col="obs_motor", flow_col="flow")
        k_mean = float(obs.sum() / max(raw_f.sum(), 1.0))
    else:
        raise ValueError(f"Unknown scaling_mode: {scaling_mode}")

    tbl["pred_motor"] = pred
    geh = np.sqrt(2.0 * (0.1 * pred - 0.1 * obs) ** 2 / np.maximum(0.1 * pred + 0.1 * obs, 1e-9))
    tbl["geh"] = geh

    r_p = float(stats.pearsonr(obs, pred).statistic) if np.ptp(pred) > 0 else np.nan
    rho = float(stats.spearmanr(obs, pred).statistic) if np.ptp(pred) > 0 else np.nan
    mae = float(np.mean(np.abs(pred - obs)))
    rmse = float(np.sqrt(np.mean((pred - obs) ** 2)))
    mape = float(np.mean(np.abs(pred - obs) / np.maximum(obs, 1.0)) * 100.0)
    med_ape = float(np.median(np.abs(pred - obs) / np.maximum(obs, 1.0)) * 100.0)
    geh5 = float(np.mean(geh < 5.0) * 100.0)
    geh10 = float(np.mean(geh < 10.0) * 100.0)
    mean_geh = float(np.mean(geh))
    sst_raw = float(np.sum((obs - obs.mean()) ** 2))
    raw_r2 = float(1.0 - np.sum((obs - pred) ** 2) / sst_raw) if sst_raw > 0 else np.nan
    lo, lp = np.log1p(obs), np.log1p(np.maximum(pred, 0.0))
    sst_log = float(np.sum((lo - lo.mean()) ** 2))
    log_r2 = float(1.0 - np.sum((lo - lp) ** 2) / sst_log) if sst_log > 0 else np.nan

    # Leeds LAD network MVKM/day
    lad_f = net_f_scaled[lad_in_mask]
    mvkm_all = float((lad_f * lad_len_m).sum() / 1e9)
    maj_mask = net_enr.loc[lad_in_mask, "road_class"].isin(["motorway", "trunk", "primary"]).to_numpy(bool)
    mvkm_maj = float((lad_f[maj_mask] * lad_len_m[maj_mask]).sum() / 1e9)

    metrics = {
        "n_sites_total": int(len(s)),
        "n_sites_matched": int(len(tbl)),
        "coverage_pos_pct": float(np.mean(raw_f > 0) * 100.0),
        "pearson_r": r_p,
        "spearman_rho": rho,
        "mae": mae,
        "rmse": rmse,
        "mape_pct": mape,
        "median_ape_pct": med_ape,
        "geh_lt_5_pct": geh5,
        "geh_lt_10_pct": geh10,
        "mean_geh": mean_geh,
        "raw_r2": raw_r2,
        "log_r2": log_r2,
        "scale_factor_k": k_mean,
        "lad_mvkm_day_all": mvkm_all,
        "lad_mvkm_day_major": mvkm_maj,
    }
    return metrics, tbl, net_f_scaled


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lad_poly = gpd.read_file(DATA_DIR / "leeds-lad.geojson").to_crs(27700).geometry.iloc[0]
    cx, cy = Transformer.from_crs(4326, 27700, always_xy=True).transform(CENTRE_LON, CENTRE_LAT)

    # Load OD & T1/T2 covariates
    zones = gpd.read_parquet(CASE_DIR / "inputs/zones.parquet")
    od_raw = pd.read_parquet(CASE_DIR / "inputs/od_synthetic.parquet")
    od_motor = od_raw[od_raw["mode"].isin(["car", "heavy"])].groupby(["o_zone", "d_zone"], as_index=False)["trips"].sum()
    pois = pd.read_parquet(CASE_DIR / "inputs/pois.parquet")
    pop = pd.read_parquet(CASE_DIR / "inputs/population.parquet")

    # Evaluation count sets
    dft_lad = pd.read_parquet(OUT_DIR / "dft_counts_leeds.parquet")
    eval_2021 = dft_lad[dft_lad["has_2021_2025"]].copy().assign(obs_motor=lambda d: d["p2021_motor"]).reset_index(drop=True)
    eval_2015 = dft_lad[dft_lad["has_2015_2025"]].copy().assign(obs_motor=lambda d: d["p2015_motor"]).reset_index(drop=True)

    bench_sites = pd.read_csv(CASE_DIR / "sites.csv")
    bench_mot = bench_sites[bench_sites["mode"].isin(["car", "heavy"])].copy()
    eval_bench = (
        bench_mot.groupby(["site_id", "x", "y", "group_id", "fold"], as_index=False)["value"]
        .sum()
        .rename(columns={"value": "obs_motor"})
    )

    net_names = [
        "N0_raw_unfiltered",
        "N1_raw_motor",
        "N2_minimal_python",
        "N3_csis_netclean",
        "N4_csis_strategic",
    ]

    all_metric_rows = []
    cat_rows = []
    site_pred_frames = []
    central_lad_gdf = None

    for n_id in net_names:
        print(f"\n--- Running assignments on {n_id} ---")
        net_gdf = gpd.read_parquet(WORK_DIR / f"{n_id}.parquet")
        net_enr = enrich_network_attributes(net_gdf, pois, pop)

        lad_in_mask = net_enr.intersects(lad_poly).to_numpy(bool)
        lad_len_m = net_enr.loc[lad_in_mask].geometry.intersection(lad_poly).length.to_numpy(float)

        cw_2021 = build_crosswalk(net_enr, eval_2021)
        cw_2015 = build_crosswalk(net_enr, eval_2015)
        cw_bench = build_crosswalk(net_enr, eval_bench)

        edges_base, od_base = build_baseline_edges_and_od(net_enr, zones, od_motor)

        # 1. AON baseline
        seg_aon, dt_aon, sc_aon, _ = run_single_class_aeq(edges_base, od_base, algorithm="all-or-nothing", max_iter=1)
        seg_aon_mod = modulate_aadtbench(net_enr, seg_aon)

        # 2. UE BFW baseline (with aadtbench 0.09 * cap_km cliff)
        seg_ue_cliff, dt_ue_c, sc_ue_c, rg_ue_c = run_single_class_aeq(
            edges_base, od_base, algorithm="bfw", max_iter=20, rgap=1e-3, use_cap_cliff=True
        )
        seg_ue_mod = modulate_aadtbench(net_enr, seg_ue_cliff)

        # 3. UE BFW nocliff (fixing ONLY the 0.09 * cap_km cliff so congestion UE stays active)
        seg_ue_nc, dt_ue_nc, sc_ue_nc, rg_ue_nc = run_single_class_aeq(
            edges_base, od_base, algorithm="bfw", max_iter=20, rgap=1e-3, use_cap_cliff=False
        )

        run_specs = [
            ("aon_physical", "AON (physical, OOF vol-scaled)", seg_aon, dt_aon, sc_aon, 0.0, "oof_vol"),
            ("aon_aadtbench", "AON (aadtbench modulated)", seg_aon_mod, dt_aon, sc_aon, 0.0, "oof_geom"),
            ("ue_bfw_physical", "UE-BFW baseline (with cap_km cliff)", seg_ue_cliff, dt_ue_c, sc_ue_c, rg_ue_c, "oof_vol"),
            ("ue_bfw_aadtbench", "UE-BFW (aadtbench modulated)", seg_ue_mod, dt_ue_c, sc_ue_c, rg_ue_c, "oof_geom"),
            ("ue_bfw_nocliff", "UE-BFW (cap_km cliff removed)", seg_ue_nc, dt_ue_nc, sc_ue_nc, rg_ue_nc, "oof_vol"),
        ]

        # On N4_csis_strategic also run Ablation A4 (netweld lanes/speeds on internal single-class OD)
        if n_id == "N4_csis_strategic":
            cap_map = net_enr.set_index("segment_id")["cap_vph"]
            spd_map = net_enr.set_index("segment_id")["speed_kmh"]
            c_arr = edges_base["segment_id"].map(cap_map).fillna(450.0).to_numpy(float)
            t_arr = edges_base["length_m"].to_numpy(float) / (edges_base["segment_id"].map(spd_map).fillna(30.0).to_numpy(float) / 3.6)
            seg_a4, dt_a4, sc_a4, rg_a4 = run_single_class_aeq(
                edges_base, od_base, algorithm="bfw", max_iter=20, rgap=1e-3, use_cap_cliff=False, custom_cap=c_arr, custom_time=t_arr
            )
            run_specs.append(("ue_welded_internal", "UE-BFW + netweld lanes & maxspeed (internal OD)", seg_a4, dt_a4, sc_a4, rg_a4, "oof_vol"))

        # 4. Improved multi-class SRN-cordon BFW UE (run on N1..N4)
        seg_imp = None
        if n_id != "N0_raw_unfiltered":
            seg_imp, dt_imp, rg_imp, n_gw = run_improved_multiclass_ue(net_enr, zones, od_raw, cx, cy)
            run_specs.append(("ue_improved_phys", "UE-Improved (physical multi-class + SRN cordon, k=1)", seg_imp, dt_imp, 0.085, rg_imp, "none"))
            run_specs.append(("ue_improved_phys_scaled", "UE-Improved (physical + 1-scalar OOF vol scale)", seg_imp, dt_imp, 0.085, rg_imp, "oof_vol"))
            run_specs.append(("ue_improved_cal", "UE-Improved + CSIS spatial OOF calibration", seg_imp, dt_imp, 0.085, rg_imp, "oof_csis_cal"))

        for var_id, var_label, seg_df, wall_s, eff_phf, rgap_val, sc_mode in run_specs:
            for eval_id, eval_df, cw_df in [
                ("leeds_lad_2021_2025", eval_2021, cw_2021),
                ("leeds_lad_2015_2025", eval_2015, cw_2015),
                ("aadtbench_leeds_v1", eval_bench, cw_bench),
            ]:
                # Skip aadtbench_leeds_v1 on oof_csis_cal if not needed? Actually it has fold & value too, runs in 5ms!
                m, tbl_out, net_scaled = evaluate_predictions(eval_df, cw_df, net_enr, seg_df, lad_in_mask, lad_len_m, scaling_mode=sc_mode)
                row = {
                    "network_id": n_id,
                    "assignment_id": var_id,
                    "assignment_label": var_label,
                    "eval_set": eval_id,
                    "wall_time_s": round(wall_s, 2),
                    "effective_phf": round(eff_phf, 5),
                    "final_rgap": rgap_val,
                    **m,
                }
                all_metric_rows.append(row)

                if eval_id == "leeds_lad_2021_2025":
                    print(
                        f"  {var_id:24s} | n={m['n_sites_matched']:3d} r={m['pearson_r']:.3f} rho={m['spearman_rho']:.3f} "
                        f"MAE={m['mae']:6.0f} MedAPE={m['median_ape_pct']:5.1f}% GEH<5={m['geh_lt_5_pct']:5.1f}% "
                        f"MVKM={m['lad_mvkm_day_all']:5.2f}M (Maj={m['lad_mvkm_day_major']:5.2f}M)"
                    )
                    # Save site-level predictions & road-category breakdown for key configurations
                    if var_id in {"aon_physical", "aon_aadtbench", "ue_bfw_physical", "ue_bfw_nocliff", "ue_improved_phys", "ue_improved_cal"}:
                        keep_s = [
                            c for c in ["site_id", "road_name", "road_type", "road_category", "x", "y", "fold", "obs_motor", "flow", "pred_motor", "geh"]
                            if c in tbl_out.columns
                        ]
                        sp = tbl_out[keep_s].copy()
                        sp["network_id"] = n_id
                        sp["assignment_id"] = var_id
                        site_pred_frames.append(sp)

                        if "road_category" in tbl_out.columns:
                            cb = (
                                tbl_out.groupby("road_category", as_index=False)
                                .agg(
                                    n_sites=("obs_motor", "count"),
                                    obs_mean=("obs_motor", "mean"),
                                    pred_mean=("pred_motor", "mean"),
                                    mae=("pred_motor", lambda s_p: float(np.mean(np.abs(s_p - tbl_out.loc[s_p.index, "obs_motor"])))),
                                    geh_lt_5_pct=("geh", lambda g: float(np.mean(g < 5.0) * 100.0)),
                                )
                            )
                            cb["network_id"] = n_id
                            cb["assignment_id"] = var_id
                            cat_rows.append(cb)

        # Build central Leeds LAD link-level AADT GeoParquet on N4_csis_strategic
        if n_id == "N4_csis_strategic" and seg_imp is not None:
            print("\n--- Building central Leeds LAD AADT GeoParquet (N4_csis_strategic) ---")
            # Calibrate motor, light (car), and heavy separately using Leeds LAD 2021-2025 sites
            s_cal = eval_2021.drop(columns=[c for c in ("motor", "car", "heavy") if c in eval_2021.columns]).copy()
            s_cal["site_id"] = s_cal["site_id"].astype(str)
            s_cal["obs_car"] = s_cal["p2021_car"]
            s_cal["obs_heavy"] = s_cal["p2021_heavy"]

            xw_imp = cw_2021.merge(seg_imp.assign(segment_id=seg_imp["segment_id"].astype(str)), on="segment_id", how="left")
            for c_fl in ("flow", "car", "heavy"):
                xw_imp[c_fl] = xw_imp[c_fl].fillna(0.0).clip(lower=0.0) * xw_imp["weight"].astype(float)
            agg_imp = xw_imp.groupby("site_id", as_index=False)[["flow", "car", "heavy"]].sum()

            prim4 = cw_2021[cw_2021["rule"] != "pair"].merge(
                net_enr[["segment_id", "road_class", "cap_vph", "speed_kmh", "log_attr"]].assign(segment_id=net_enr["segment_id"].astype(str)),
                on="segment_id", how="left",
            )
            cap2w_4 = (
                cw_2021.merge(net_enr[["segment_id", "cap_vph"]].assign(segment_id=net_enr["segment_id"].astype(str)), on="segment_id", how="left")
                .groupby("site_id")["cap_vph"]
                .sum()
            )
            tbl4 = s_cal.merge(agg_imp, on="site_id", how="inner")
            tbl4["road_class"] = tbl4["site_id"].map(prim4.set_index("site_id")["road_class"]).fillna("minor")
            tbl4["cap_2w"] = tbl4["site_id"].map(cap2w_4).fillna(900.0)
            tbl4["speed_kmh"] = tbl4["site_id"].map(prim4.set_index("site_id")["speed_kmh"]).fillna(30.0)
            tbl4["log_attr"] = tbl4["site_id"].map(prim4.set_index("site_id")["log_attr"]).fillna(0.0)

            oof_mot, net_mot_cal = fit_oof_and_network_calibration(tbl4, net_enr, seg_imp, target_col="obs_motor", flow_col="flow")
            oof_car, net_car_cal = fit_oof_and_network_calibration(tbl4, net_enr, seg_imp, target_col="obs_car", flow_col="car")
            oof_hvy, net_hvy_cal = fit_oof_and_network_calibration(tbl4, net_enr, seg_imp, target_col="obs_heavy", flow_col="heavy")

            # Reconcile car + heavy shares proportionally so aadt_car + aadt_heavy == aadt_motor identically
            sum_parts = np.maximum(net_car_cal + net_hvy_cal, 1e-6)
            net_car_rec = net_mot_cal * (net_car_cal / sum_parts)
            net_hvy_rec = net_mot_cal * (net_hvy_cal / sum_parts)

            print(
                f"  Sub-mode OOF Pearson r -> Light (car+LGV): {stats.pearsonr(tbl4['obs_car'], oof_car).statistic:.3f} | "
                f"Heavy (HGV+bus): {stats.pearsonr(tbl4['obs_heavy'], oof_hvy).statistic:.3f}"
            )

            lad_links = net_enr.loc[lad_in_mask].copy().reset_index(drop=True)
            lad_links["lad_length_m"] = np.round(lad_len_m, 1)
            phys_merge = (
                lad_links[["segment_id"]]
                .assign(segment_id=lad_links["segment_id"].astype(str))
                .merge(seg_imp.assign(segment_id=seg_imp["segment_id"].astype(str)), on="segment_id", how="left")
                .fillna(0.0)
            )
            lad_links["aadt_motor"] = np.round(net_mot_cal[lad_in_mask], 1)
            lad_links["aadt_car"] = np.round(net_car_rec[lad_in_mask], 1)
            lad_links["aadt_heavy"] = np.round(net_hvy_rec[lad_in_mask], 1)
            lad_links["aadt_motor_phys"] = np.round(phys_merge["flow"].to_numpy(float), 1)
            lad_links["aadt_car_phys"] = np.round(phys_merge["car"].to_numpy(float), 1)
            lad_links["aadt_heavy_phys"] = np.round(phys_merge["heavy"].to_numpy(float), 1)
            lad_links["vc_ratio_peak"] = np.round(phys_merge["vc_ratio"].to_numpy(float), 3)
            lad_links["cap_vph"] = np.round(lad_links["cap_vph"].to_numpy(float), 1)
            lad_links["speed_kmh"] = np.round(lad_links["speed_kmh"].to_numpy(float), 1)

            # Simplify geometry to 0.5 m tolerance & 0.1 m coordinate precision to keep committed parquet compact (<6 MB)
            lad_links["geometry"] = shapely.set_precision(shapely.simplify(lad_links.geometry.values, 0.5), 0.1)
            keep_lad_cols = [
                "segment_id", "road_class", "highway", "name", "ref", "dir_car",
                "length_m", "lad_length_m", "dir_lanes", "cap_vph", "speed_kmh",
                "aadt_motor", "aadt_car", "aadt_heavy",
                "aadt_motor_phys", "aadt_car_phys", "aadt_heavy_phys", "vc_ratio_peak",
                "geometry",
            ]
            keep_lad_cols = [c for c in keep_lad_cols if c in lad_links.columns]
            central_lad_gdf = lad_links[keep_lad_cols]

    metrics_df = pd.DataFrame(all_metric_rows)
    metrics_df.to_parquet(OUT_DIR / "validation_metrics.parquet", index=False)

    cat_df = pd.concat(cat_rows, ignore_index=True)
    cat_df.to_parquet(OUT_DIR / "category_breakdown.parquet", index=False)

    site_preds_df = pd.concat(site_pred_frames, ignore_index=True)
    site_preds_df.to_parquet(OUT_DIR / "site_predictions.parquet", index=False)

    if central_lad_gdf is not None:
        central_lad_gdf.to_parquet(OUT_DIR / "leeds_lad_aadt.parquet", index=False, compression="zstd")

    print("\nSaved all validation and central estimate parquets to data/processed/.")


if __name__ == "__main__":
    main()

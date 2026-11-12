#!/usr/bin/env python3
"""01_prepare_counts.py — Prepare DfT AADF motor traffic counts and VKM reference for Leeds LAD.

Reads `data/leeds-lad.geojson` and `data/dft/*.parquet` (normalised by `countsget`),
aggregates motor classes (`motor`, `car` [light motor], `heavy`), constructs spatial
block cross-validation folds (3 km grid, 5 folds, seed=42), and computes the
DfT-implied Leeds LAD motor-vehicle-kilometres/day (MVKM) reference totals.

Outputs:
  - data/processed/dft_counts_leeds.parquet
  - data/processed/dft_vkm_reference.parquet
"""
from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
DFT_DIR = DATA_DIR / "dft"
OUT_DIR = DATA_DIR / "processed"


def assign_spatial_folds(df: pd.DataFrame, block_m: float = 3000.0, n_folds: int = 5, seed: int = 42) -> pd.DataFrame:
    """Assign deterministic 3 km spatial block folds (matching aadtbench spatial CV)."""
    out = df.copy()
    gx = np.floor(out["x"].to_numpy(float) / block_m).astype(int)
    gy = np.floor(out["y"].to_numpy(float) / block_m).astype(int)
    out["group_id"] = "g" + pd.Series(gx, index=out.index).astype(str) + "_" + pd.Series(gy, index=out.index).astype(str)
    groups = np.array(sorted(out["group_id"].unique()))
    order = np.random.default_rng(seed).permutation(len(groups))
    fold_map = {groups[g]: int(rank % n_folds) for rank, g in enumerate(order)}
    out["fold"] = out["group_id"].map(fold_map).astype(int)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load Leeds LAD boundary in EPSG:27700
    lad_gdf = gpd.read_file(DATA_DIR / "leeds-lad.geojson").to_crs(27700)
    lad_poly = lad_gdf.geometry.iloc[0]
    lad_area_km2 = float(lad_poly.area / 1e6)

    # 2. Load DfT sites and unpack JSON metadata
    sites = pd.read_parquet(DFT_DIR / "sites.parquet")
    meta = pd.json_normalize(sites["metadata"].apply(json.loads))
    sites = pd.concat([sites.drop(columns=["metadata"]), meta], axis=1)
    sites["link_length_km"] = pd.to_numeric(sites.get("link_length_km"), errors="coerce")
    sites["link_length_miles"] = pd.to_numeric(sites.get("link_length_miles"), errors="coerce")

    gsites = gpd.GeoDataFrame(
        sites,
        geometry=gpd.points_from_xy(sites["lon"], sites["lat"]),
        crs=4326,
    ).to_crs(27700)
    gsites["x"] = gsites.geometry.x
    gsites["y"] = gsites.geometry.y
    gsites["in_leeds_lad"] = gsites.within(lad_poly)

    # 3. Load channels & annual mean_daily observations
    channels = pd.read_parquet(DFT_DIR / "channels.parquet")
    obs = pd.read_parquet(DFT_DIR / "observations.parquet")
    oc = obs.merge(channels[["channel_id", "site_id", "class", "measure"]], on="channel_id", how="inner")
    oc = oc[oc["measure"] == "mean_daily"].copy()
    oc["year"] = pd.to_datetime(oc["period_start"]).dt.year

    # Classify vehicle groups:
    # motor = cars/taxis + two-wheeled motor + LGVs + buses/coaches + rigid/articulated HGVs
    is_motor = oc["class"].str.startswith("motor.")
    oc_motor = oc[is_motor].copy()
    oc_motor["veh_group"] = np.where(
        oc_motor["class"].str.startswith("motor.heavy"),
        "heavy",
        np.where(oc_motor["class"] == "motor.car", "car_only", "other_light"),
    )

    annual = (
        oc_motor.groupby(["site_id", "year", "veh_group"])["count"]
        .sum()
        .unstack(fill_value=0.0)
        .reset_index()
    )
    for col in ("car_only", "other_light", "heavy"):
        if col not in annual.columns:
            annual[col] = 0.0
    # `car` denotes light motor vehicles (cars/taxis + LGVs + two-wheelers); `car_only` is cars/taxis alone
    annual["car"] = annual["car_only"] + annual["other_light"]
    annual["motor"] = annual["car"] + annual["heavy"]

    # Cross-check 2021-2025 mean against summary_dft.parquet
    summary_dft = pd.read_parquet(DFT_DIR / "summary_dft.parquet")
    sum_piv = summary_dft.pivot(index="site_id", columns="group", values="mean_daily").fillna(0.0)
    sum_piv["summary_motor"] = sum_piv.get("car", 0.0) + sum_piv.get("heavy", 0.0) + sum_piv.get("other", 0.0)

    def period_summary(df_ann: pd.DataFrame, min_year: int, prefix: str) -> pd.DataFrame:
        sub = df_ann[df_ann["year"] >= min_year]
        agg = (
            sub.groupby("site_id")
            .agg(
                **{
                    f"{prefix}_motor": ("motor", "mean"),
                    f"{prefix}_car": ("car", "mean"),
                    f"{prefix}_car_only": ("car_only", "mean"),
                    f"{prefix}_heavy": ("heavy", "mean"),
                    f"{prefix}_n_years": ("year", "nunique"),
                    f"{prefix}_last_year": ("year", "max"),
                }
            )
            .reset_index()
        )
        return agg

    p_2021 = period_summary(annual, 2021, "p2021")
    p_2015 = period_summary(annual, 2015, "p2015")
    p_all = period_summary(annual, 2000, "pall")

    # Verify exact agreement with summary_dft.parquet
    chk = p_2021.merge(sum_piv[["summary_motor"]].reset_index(), on="site_id", how="inner")
    max_diff = float(np.max(np.abs(chk["p2021_motor"] - chk["summary_motor"])))
    assert max_diff < 1e-6, f"2021-2025 motor AADF mismatch vs summary_dft: {max_diff}"

    # 4. Merge site attributes andperiod summaries (restrict to sites inside Leeds LAD + keep flag)
    keep_cols = [
        "site_id", "source_id", "name", "lon", "lat", "x", "y", "in_leeds_lad",
        "road_name", "road_type", "road_category", "local_authority_name",
        "link_length_km", "start_junction_road_name", "end_junction_road_name",
    ]
    keep_cols = [c for c in keep_cols if c in gsites.columns]
    merged = (
        pd.DataFrame(gsites[keep_cols])
        .merge(p_2021, on="site_id", how="left")
        .merge(p_2015, on="site_id", how="left")
        .merge(p_all, on="site_id", how="left")
    )

    # Assign 5-fold 3km spatial block CV separately on the primary Leeds LAD 2021-2025 set
    # and on the extended 2015-2025 Leeds LAD set so every site has a deterministic block & fold
    lad_sites = merged[merged["in_leeds_lad"] & merged["pall_motor"].notna()].copy().reset_index(drop=True)
    lad_sites = assign_spatial_folds(lad_sites, block_m=3000.0, n_folds=5, seed=42)

    # Primary columns (`motor`, `car`, `heavy`) use 2021-2025 post-COVID mean where available
    lad_sites["has_2021_2025"] = lad_sites["p2021_motor"].notna()
    lad_sites["has_2015_2025"] = lad_sites["p2015_motor"].notna()
    lad_sites["motor"] = lad_sites["p2021_motor"]
    lad_sites["car"] = lad_sites["p2021_car"]
    lad_sites["heavy"] = lad_sites["p2021_heavy"]

    lad_sites.to_parquet(OUT_DIR / "dft_counts_leeds.parquet", index=False)

    # 5. Compute DfT-implied network motor-vehicle-km/day (MVKM) reference for Leeds LAD
    # Pre- vs post-COVID adjustment factor on overlapping Leeds LAD Major road sites
    pre_2015_19 = (
        annual[(annual["year"] >= 2015) & (annual["year"] <= 2019)]
        .groupby("site_id")[["motor", "car", "heavy"]]
        .mean()
        .rename(columns=lambda c: f"pre_{c}")
        .reset_index()
    )
    overlap = lad_sites[lad_sites["road_type"] == "Major"].merge(pre_2015_19, on="site_id", how="inner")
    overlap = overlap[overlap["p2021_motor"].notna()]
    post_pre_ratio = float(overlap["p2021_motor"].sum() / overlap["pre_motor"].sum())

    maj = lad_sites[(lad_sites["road_type"] == "Major") & lad_sites["link_length_km"].notna()].copy()
    # Use 2021-2025 AADF directly where observed (150 of 211 in-boundary Major links),
    # and post-COVID ratio-adjusted historical AADF on the remaining links so all 394 km of Leeds MRN are covered
    maj["adj_motor"] = np.where(maj["p2021_motor"].notna(), maj["p2021_motor"], maj["pall_motor"] * post_pre_ratio)
    maj["adj_car"] = np.where(maj["p2021_car"].notna(), maj["p2021_car"], maj["pall_car"] * post_pre_ratio)
    maj["adj_heavy"] = np.where(maj["p2021_heavy"].notna(), maj["p2021_heavy"], maj["pall_heavy"] * post_pre_ratio)
    maj["mvkm_day_raw_latest"] = maj["pall_motor"] * maj["link_length_km"] / 1e6
    maj["mvkm_day_postcovid"] = maj["adj_motor"] * maj["link_length_km"] / 1e6
    maj["car_vkm_day_postcovid"] = maj["adj_car"] * maj["link_length_km"] / 1e6
    maj["heavy_vkm_day_postcovid"] = maj["adj_heavy"] * maj["link_length_km"] / 1e6

    by_cat = (
        maj.groupby("road_category", as_index=False)
        .agg(
            n_links=("site_id", "count"),
            n_postcovid_observed=("has_2021_2025", "sum"),
            road_km=("link_length_km", "sum"),
            mean_aadf_postcovid=("adj_motor", "mean"),
            mvkm_day_postcovid=("mvkm_day_postcovid", "sum"),
            mvkm_day_raw_latest=("mvkm_day_raw_latest", "sum"),
            car_vkm_day=("car_vkm_day_postcovid", "sum"),
            heavy_vkm_day=("heavy_vkm_day_postcovid", "sum"),
        )
    )
    tot_row = pd.DataFrame(
        [
            {
                "road_category": "ALL_MAJOR",
                "n_links": int(by_cat["n_links"].sum()),
                "n_postcovid_observed": int(by_cat["n_postcovid_observed"].sum()),
                "road_km": float(by_cat["road_km"].sum()),
                "mean_aadf_postcovid": float(
                    (maj["adj_motor"] * maj["link_length_km"]).sum() / maj["link_length_km"].sum()
                ),
                "mvkm_day_postcovid": float(by_cat["mvkm_day_postcovid"].sum()),
                "mvkm_day_raw_latest": float(by_cat["mvkm_day_raw_latest"].sum()),
                "car_vkm_day": float(by_cat["car_vkm_day"].sum()),
                "heavy_vkm_day": float(by_cat["heavy_vkm_day"].sum()),
            }
        ]
    )
    # DfT published LA-level total traffic reference (TRA8904 Leeds all motor vehicles ~4.23 bn vehicle-miles/yr = 18.65 M veh-km/day)
    # Minor roads share in urban Yorkshire metropolitan districts ~34-36% of total MVKM
    maj_mvkm = float(tot_row["mvkm_day_postcovid"].iloc[0])
    est_all_roads = pd.DataFrame(
        [
            {
                "road_category": "IMPLIED_ALL_ROADS_DFT_SHARE",
                "n_links": int(len(lad_sites)),
                "n_postcovid_observed": int(lad_sites["has_2021_2025"].sum()),
                "road_km": np.nan,
                "mean_aadf_postcovid": np.nan,
                "mvkm_day_postcovid": maj_mvkm / 0.65,
                "mvkm_day_raw_latest": float(tot_row["mvkm_day_raw_latest"].iloc[0]) / 0.65,
                "car_vkm_day": np.nan,
                "heavy_vkm_day": np.nan,
            }
        ]
    )
    vkm_ref = pd.concat([by_cat, tot_row, est_all_roads], ignore_index=True)
    vkm_ref["lad_area_km2"] = lad_area_km2
    vkm_ref["post_pre_covid_ratio"] = post_pre_ratio
    vkm_ref.to_parquet(OUT_DIR / "dft_vkm_reference.parquet", index=False)

    print(f"Leeds LAD area: {lad_area_km2:.2f} km^2")
    print(
        f"Leeds LAD DfT sites: {len(lad_sites)} all-years | "
        f"{int(lad_sites['has_2015_2025'].sum())} (2015-2025) | "
        f"{int(lad_sites['has_2021_2025'].sum())} (2021-2025 primary)"
    )
    print(f"Summary_dft exact check max diff: {max_diff:.3e} | Post/Pre-COVID Major ratio: {post_pre_ratio:.4f}")
    print("\nDfT Major Road Network VKM Reference (Leeds LAD):")
    print(vkm_ref[["road_category", "n_links", "n_postcovid_observed", "road_km", "mvkm_day_postcovid", "mvkm_day_raw_latest"]].to_string(index=False))


if __name__ == "__main__":
    main()

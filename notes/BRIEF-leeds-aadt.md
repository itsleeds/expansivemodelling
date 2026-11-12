# Brief: Leeds motor traffic AADT with AequilibraE — expansivemodelling

You are working in the repository `~/github/itsleeds/expansivemodelling` (your cwd),
a new public repo. Your job is to build a central estimate of motor traffic AADT
across the Leeds local authority district network with AequilibraE, validate it
against DfT counts, quantify the benefit of network cleaning, and write the
results up as a camera-ready Quarto paper with reproducible code chunks.

## Inputs (all present)

- `data/leeds-lad.geojson` — Leeds LAD boundary (ONS, LAD24CD E08000035).
- `data/dft/` — DfT AADF traffic counts normalised by `countsget`
  (atf-lab/countsget): `sites_dft.parquet`, `channels_dft.parquet`,
  `observations_dft.parquet`, `daily_dft.parquet`, `summary_dft.parquet`.
  AADF rows have `measure = mean_daily`; the bin is the year. Aggregate motor
  classes (cars + two-wheeled motor + buses/coaches + light/heavy goods) into
  `motor` (and keep `car` and `heavy` separately if useful).
- Raw OSM network and count sites from the aadtbench benchmark case
  `~/github/robinlovelace/aadtbench/cases/leeds-v1/`: `network.parquet`
  (raw OSM links, ~230k links, no cleaning), `source.osm.pbf`, `sites.csv`,
  `crosswalk.csv`, and `inputs/` (od_synthetic.parquet, zones.parquet,
  population.parquet, pois.parquet). Read
  `~/github/robinlovelace/aadtbench/BENCHMARK.md` and `adapters/_base.py` for
  the network schema and `road_class` definitions.

## Starting settings (from aadtbench, which worked well)

Start from `~/github/robinlovelace/aadtbench/adapters/aequilibrae.py`:

- Variant `aon`: all-or-nothing on free-flow time, distance="time".
- Variant `ue_bfw`: BFW user equilibrium, max_iter 20, rgap 1e-3, BPR VDF with
  alpha=0.15, beta=4.0; demand scaled by peak hour factor (daily trips x 0.1,
  flows x 10 back) with the capacity-based demand scale check (cap in
  vph from `CLASS_CAPACITY_VPH` mapped by `road_class`, default 300 vph).
- OD demand: `inputs/od_synthetic.parquet` (T2 synthetic OD tier).
- You may improve on these settings, but report the aadtbench settings as the
  baseline and keep the comparison: baseline (aon), congestion-aware (ue_bfw),
  plus your improved configuration.

## Network cleaning study (important)

The aim includes showing the benefits of network cleaning using open tools from
`~/github`, especially those built in the csis project. Survey and use what
fits (read their READMEs first; they are read-only for you):

- `~/github/robinlovelace/csis` and `~/github/robinlovelace/criticalissues`
  (Rust crates + pipeline; look for network cleaning, simplification,
  junction handling).
- `~/github/atf-lab/osmactive`, `~/github/atf-lab/netfit`,
  `~/github/robinlovelace/netassign`, `~/github/atfutures/netweld`,
  `~/github/robinlovelace/osmget`.

Build a comparison of at least: (1) raw OSM as-is, (2) a minimally cleaned
network (motor-accessible filter, largest component, dangling-link trimming),
(3) a further-cleaned network using a csis/osmactive-style procedure if
available. Report per variant: link count, network km, component count,
largest-component share, and the validation metrics below. State clearly which
tool produced which cleaning.

## Validation

- Join modelled link AADT to DfT AADF count sites (spatial match of count
  points to network links; use `sites.csv` crosswalk from the benchmark case
  as a starting point and document the matching rule).
- Report: number of matched sites, Pearson/Spearman rho, MAE, MAPE, GEH
  statistics (share of sites with GEH < 5), and the network-level total
  motor-vehicle km or total daily flow vs the DfT-implied total.
- Include the uncritical baseline: raw network + aon.

## Deliverables

1. `paper/index.qmd` — the paper. Camera-ready: title, abstract, methods
   (inputs, settings, cleaning procedures with tool attributions), results
   (tables + figures), discussion, references. All chunks reproducible:
   - Compute chunks read committed caches under `data/processed/`; heavy
     generation scripts live in `scripts/` and are documented with their exact
     invocation in an appendix chunk (eval=false) and in the README.
   - Total `quarto render` wall time under ~10 minutes with caches committed.
   - Figures under `paper/figs/` (committed), saved by the scripts.
2. `scripts/` — the actual build scripts (Python; create a repo venv at
   `.venv` and `pip install -r requirements.txt`).
3. `data/processed/` — committed parquet caches small enough for git (<50 MB
   total); anything larger stays uncommitted and gitignored.
4. Update the repo README: how to reproduce end to end.
5. Commit your work to a new branch `agy/leeds-aadt`, pushing to origin
   (`git push -u origin agy/leeds-aadt`) as you reach milestones, and finish by
   opening a PR to main (`gh pr create`) with a summary of what was done and
   what did NOT work.

## Constraints

- Loose lease: you may run long jobs, commit and push to your branch freely;
  do NOT merge, do NOT push to main, do NOT modify files outside this repo.
- Do not re-download large datasets that already exist; use the inputs above.
- OpenStreetMap-derived outputs: ODbL attribution; DfT: OGL v3 — include in
  the paper.
- If a step fails after honest attempts, document it in the PR rather than
  fabricating results. Numbers in the paper must come from code you ran.

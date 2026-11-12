# expansivemodelling

Reproducible motor traffic modelling for Leeds with open tools, benchmarked
against DfT traffic counts. The starting point is the AequilibraE traffic
assignment configuration that performed well in the
[aadtbench](https://github.com/Robinlovelace/aadtbench) benchmark
(`leeds-v1` case, `aon` and `ue_bfw` variants).

## Aims

1. Generate a central estimate of motor traffic AADT across the Leeds local
   authority district network, using cleaned OpenStreetMap data as input.
2. Quantify the benefits of network cleaning (topology simplification,
   junction consolidation) using open tools from `~/github`, especially those
   developed in the [csis](https://github.com/Robinlovelace/csis) project.
3. Compare open tools (AequilibraE, cityseer, sDNA, osmactive/netweld) with
   closed tools (SATURN and similar proprietary packages) in a documented
   feature table.

## Repository layout

- `paper/index.qmd` — the main paper, with reproducible code chunks that
  produce camera-ready figures and tables.
- `features_table.qmd` — feature comparison of open vs closed traffic
  modelling tools.
- `data/` — committed inputs: Leeds LAD boundary and DfT traffic counts
  (obtained with [`countsget`](https://github.com/atf-lab/countsget)).
- `.github/workflows/site.yml` — renders the Quarto site to GitHub Pages on
  push (same pattern as aadtbench).

## Data inputs

| Input | Source | Path |
|---|---|---|
| Leeds LAD boundary | ONS Open Geography Portal (OGL v3) | `data/leeds-lad.geojson` |
| DfT AADF traffic counts | DfT Road Traffic Statistics, via countsget (OGL v3) | `data/dft/` |
| Raw OSM network for Leeds | aadtbench `leeds-v1` case release | `~/github/robinlovelace/aadtbench/cases/leeds-v1/` on the build machine |

## Licences

- DfT AADF and ONS boundaries: Open Government Licence v3.0.
- OpenStreetMap: ODbL 1.0.

# Brief: Feature comparison — SATURN and closed tools vs open tools for motor traffic modelling

You are researching for the itsleeds/expansivemodelling project (repo in your cwd).
Use the web for facts about the tools; cite product/paper pages in a `## Sources` section.

## Task
Produce a rigorous feature comparison of proprietary strategic transport models
against open tools, focused on MOTOR traffic modelling (car + HGV assignment):

Closed tools: SATURN (AtkinsRéalis), PTV VISUM/Visum, PTV Vistro, EMME (Bentley),
TransCAD, CUBE (Bentley), OmniTRANS (DAT.Mobility), Aimsum/Simulation alternatives.
Open tools: AequilibraE, MATSim, OpenTRIPS/OpenSUM, FourStepModel implementations,
pandana/cityseer/sDNA (network analysis), OsRM/Valhalla routing engines as assignment bases,
netassign/netweld/osmactive (AT Futures tools, in ~/github — read their READMEs for capability notes).

## Comparison dimensions (rows of the final table)
1. Assignment methods (all-or-nothing, user equilibrium / BPR variants, stochastic, dynamic)
2. Matrix estimation / OD synthesis from counts
3. Network build & cleaning from OSM (junction handling, topology repair)
4. Zoning & demand model components (trip generation/distribution)
5. Scripting / reproducibility (CLI, Python/R API, CI renderability)
6. Licensing & cost (closed licence cost basis vs open licence terms)
7. Validation workflows (screenline/flow comparison, quality metrics)
8. Multi-resolution (macro/meso/micro coupling)
9. Anything SATURN uniquely does that open tools cannot replicate, and vice versa

## Deliverable
Write a self-contained markdown file at `notes/tool-landscape-agy.md` containing:
- a markdown table: rows = features, columns = SATURN / other closed / AequilibraE / other open,
  cells = yes/no/partial + 1-line evidence
- short prose per closed tool on what it is, its licence model, and where it excels
- what open tools are missing to fully replace SATURN-type workflows
- `## Sources` with URLs
Do NOT modify any other files in the repo. Commit nothing; leave the file on disk.

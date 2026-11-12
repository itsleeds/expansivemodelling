# Feature comparison: closed vs open tools for motor traffic modelling (car + HGV assignment)

Report prepared for the `itsleeds/expansivemodelling` features page. Scope: strategic and local
**motor traffic assignment** (car and HGV flows). All facts sourced from vendor documentation,
product pages and software repositories listed under *Sources* (retrieved November 2026);
pricing and feature claims are quoted from vendor pages and should be treated as
vendor statements. Where sources were dated, ambiguous or unavailable, this is marked
explicitly rather than guessed.

---

## 1. Comparison table

Legend: ✔ = capability present; ◐ = partial/limited; ✘ = not available; "evid." = evidence/source basis.

| # | Feature dimension | SATURN (AtkinsRéalis) | Other closed tools (VISUM / EMME / TransCAD / CUBE / OmniTRANS) | AequilibraE | Other open tools (MATSim, cityseer, sDNA, pandana, osmactive) |
|---|---|---|---|---|---|
| 1 | **Assignment: all-or-nothing** | ✔ AoN available among assignment options (SATURN help, §1.1) | ✔ All five offer AoN / incremental loading as standard (PTV FAQ; Caliper docs; Bentley CUBE docs) | ✔ AoN is the base step of every algorithm (`TrafficAssignment`, algorithm list includes `all-or-nothing`) | ◐ MATSim loads demand agent-by-agent (queue sim), not classical AoN; sDNA offers betweenness "flow models" (BtA) as a heuristic flow proxy; cityseer has origin–destination betweenness; pandana ✘ (accessibility only) |
| 2 | **Assignment: user equilibrium with BPR-style VDFs** | ✔ Wardrop UE via Frank–Wolfe family, plus Origin-Based Assignment (Bar-Gera) for high-convergence MUE (SATURN help §1.1, §21) | ✔ VISUM: Equilibrium/LUCE/BFW, DUE, ICA; EMME: FW + SOLA; TransCAD: multi-threaded high-convergence UE; CUBE Voyager HIGHWAY: UE with congestion feedback; OmniTRANS: static allocations (PTV FAQ; Bentley OpenPaths docs; Caliper docs; Bentley KB) | ✔ Multi-class UE with FW / CFW / BFW / MSA; VDFs: BPR, BPR2, CONICAL, INRETS, AKCELIK; per-class PCE, VoT and fixed costs (AequilibraE docs, *Traffic Assignment Procedure*) | ✘ None of these provide a classical static UE engine; MATSim reaches a dynamically consistent state through its iterative co-evolutionary loop, which is conceptually different from link-based Wardrop UE |
| 3 | **Stochastic / dynamic assignment** | ✔ Burrell stochastic (SUE) option; quasi-dynamic junction simulation of queue growth/decay (SATURN help §1.1) | ✔ VISUM: DUE (with spillback), DSA, simulation-based SBA mesoscopic assignment; OmniTRANS: STAQ semi-dynamic and STREAMLINE dynamic modules; EMME/CUBE/TransCAD: dynamic capabilities via DYNAMEQ / TransModeler siblings (PTV FAQ; Goudappel OmniTRANS module pages) | ✘ Static only (no stochastic or dynamic assignment in documented releases) | ◐ MATSim is fully dynamic (second-by-second queue-based mobsim) but is an agent framework, not a DTA-equilibrium module; cityseer's pedestrian-flow betweenness is stochastic-ish but pedestrian-scale |
| 4 | **Matrix estimation / OD synthesis from counts** | ✔ SATME2 matrix estimation against link counts, integrated with assignment loop (Pijα factor method) (SATURN help; UGM materials) | ✔ VISUM: TFlowFuzzy + least-squares ODME (static and dynamic); EMME: generalised-minimisation ODME; TransCAD: flexible OD estimation from counts (Nielsen method also in TransModeler); CUBE: CUBE Analyst; OmniTRANS: dedicated Matrix module (PTV FAQ; Bentley OpenPaths help; Caliper docs; Goudappel) | ◐ IPF (Iterative Proportional Fitting), synthetic gravity model calibration/apply, and a matrix-estimation procedure for matching modelled to observed flows are provided (QAequilibraE menus; docs) — but documented ODME maturity/trace record is far shorter than the commercial packages | ✘ Mostly absent; MATSim synthesises demand from population/survey data rather than estimating matrices from counts; sDNA's Two-Phase Betweenness is a generation–distribution flow heuristic, not ODME |
| 5 | **Network build from OpenStreetMap (junction handling, topology repair)** | ✘/◐ SATURN networks are traditionally hand-coded (SATCoder aids junction coding); no documented OSM importer found in public help pages (uncertainty flagged) | ◐ VISUM: documented OSM import workflow ("From OSM via PTV Visum to PTV Vissim" example); TransCAD ships routable HERE data as alternative; EMME/CUBE: import via generic GIS tools — OSM-specific junction repair not prominently documented | ✔ `create_from_osm()` from place name/bounding box/polygon; parameterisable field mapping from OSM tags (lanes split/divide, maxspeed), mode filters, unknown-tag handling, `clean` option; import from GMNS and link layers (AequilibraE docs) | ◐ cityseer includes graph-cleaning methods to reduce topological distortions; osmactive (R) extracts drive and cycling networks from OSM with infrastructure classification; MATSim needs external OSM→network conversion tools; sDNA/pandana consume standard GIS data (pandana typically fed by osmnx) |
| 6 | **Zoning & demand model components (trip generation, distribution)** | ✔ Trip matrix demand model: distribution, modal split, convergent variable-demand loops coupled to assignment; growth factoring in matrix manipulation suite (SATURN help §1.1) | ✔ Extensive: VISUM 4-step, EVA/VISEM tour-based, nested logit, ABM interfaces; CUBE GENERATION/DISTRIBUTION/FRATAR programs + CUBE Land; TransCAD "most comprehensive" demand suite incl. activity models; OmniTRANS incl. Octavius person-level demand; EMME demand adjustment tools (PTV examples page; Bentley KB; Caliper brochure; Goudappel) | ◐ IPF for balancing P/A vectors, calibrate/apply synthetic gravity model (power/negative exponential), trip length distribution reporting; no full trip-generation framework (QAequilibraE docs) | ✘/◐ MATSim takes externally generated initial demand and then optimises schedules; cityseer/sDNA/pandana do not do demand modelling; osmactive is infrastructure extraction only |
| 7 | **Scripting APIs & reproducibility (CLI, headless/CI)** | ◐ Batch operation and MSHELL manipulation shell exist (UGM presentations), but no documented Python/R API; Windows GUI-centric suite | ✔ VISUM: COM + Python automation; EMME: comprehensive Python (Modeller/Desktop/Network/Matrix APIs, Notebook); CUBE: Voyager scripting language + CubePy Python wrapper; TransCAD: GISDK Caliper Script; OmniTRANS: Python integration (vendor docs) | ✔ Pure Python/Cython package (pip-installable), native CLI added in v1.7, YAML-configured runs, OMX matrix interchange, QGIS GUI, runs headless on Linux/Windows — strong fit for CI and scripted pipelines (GitHub releases; docs) | ✔ MATSim (Java, Maven/Gradle, headless CLI), sDNA (CLI + Python API), cityseer & pandana (Python), osmactive (R) — all CI-friendly by design |
| 8 | **Licensing & cost model** | Commercial; distributed solely by Atkins(Realis) for co-owners (Van Vliet, University of Leeds, Atkins); licence fee + optional annual support; prices not published (contact-only). Windows-only (SATURN help §1.4; Wikipedia) | Commercial: VISUM from ~€9,000, sized by zones + modules (CIVITAS listing), floating/cloud via PTV Hub; OpenPaths EMME+CUBE bundles from ~$9,951 (Advanced) / ~$12,877 (Ultimate) per year (Virtuoso); TransCAD US$15,000 standard or US$6,000–7,200/yr (published price list); OmniTRANS: contact-only, 60-day 25-zone demo (vendor pages) | ✔ Free and open source (libre licence in repo — exact licence text should be checked before redistribution; repo lists "Other/NOASSERTION") | ✔ Free/open: cityseer AGPL-3.0, sDNA+ AGPL-3.0 (sDNA Open GPL-3.0), MATSim open source (GPL family — verify per module), pandana open source, osmactive open source |
| 9 | **Validation workflows (screenlines, GEH, quality metrics)** | ✔ Built around UK practice: GEH reporting of assigned vs observed flows, screenline/cordon comparisons, convergence diagnostics, SatView visualisation (SATURN help; SATURN project report literature) | ✔ TransCAD has a dedicated screenline-analysis tool and assignment-difference utilities; VISUM path analysis vs observed travel times, flows comparison; EMME/CUBE reporting toolboxes; OmniTRANS calibration workflows (Caliper docs; PTV examples; vendor pages) | ◐ Assignment produces relative-gap diagnostics, select-link analysis, flow comparisons; no packaged GEH/screenline calculator — straightforward to script in Python but not out-of-the-box | ◐ MATSim analyses/KPIs framework (counts-based calibration possible via scoring/counts modules); sDNA validates against observed flows for its heuristics; no formal GEH/screenline tooling in cityseer/pandana/osmactive |

---

## 2. The closed tools in brief

### SATURN (AtkinsRéalis)
SATURN (Simulation and Assignment of Traffic to Urban Road Networks) was developed at the
Institute for Transport Studies, University of Leeds, and has been distributed by Atkins since
1982; it remains co-owned by Van Vliet, the University of Leeds and Atkins(Realis). It is the
de-facto UK-standard congested-network assignment package, used for National Highways regional
transport models and TfL highway assignment models per User Group Meeting materials. Its
distinguishing feature is a **junction-level simulation layer built on cyclic flow profiles**,
with capacities informed by ARCADY/PICADY/OSCADY-style junction models, explicit modelling of
queue growth/decay, blocking back and flow metering between junctions. Assignment options span
AoN, Wardrop UE (Frank–Wolfe), Burrell stochastic assignment, and post-10.9 Origin-Based
Assignment for tightly converged multi-user-class assignment (cars, taxis, HGVs etc. as separate
classes). It also provides SATME2 matrix estimation, variable-demand loops (distribution, modal
split) and GPU acceleration (SatGPU). Licence fees are not published; support is sold separately.
It excels at congested urban assignment where junction interactions matter, and at UK-standard
strategic models.

### PTV VISUM
Multimodal strategic planning software from PTV (Karlsruhe; PTV Group under Porsche SE).
Static assignment: equilibrium (Wardrop), LUCE, bi-conjugate Frank–Wolfe, incremental, stochastic,
plus ICA (intersection-capacity analysis) and TRIBUT toll-sensitive assignment; dynamic: DUE with
spillback, DSA with time slices, and simulation-based assignment (SBA), a mesoscopic
vehicle-based iterative method. ODME via TFlowFuzzy and least squares (static and dynamic).
Demand side covers classic 4-step through tour-based (VISEM) and interfaces to ActivitySim.
Python/COM automation; OSM import documented. Licences sized by zone count and modules,
starting around €9,000 (CIVITAS listing), with floating/cloud options through PTV Hub.

### EMME (Bentley OpenPaths EMME)
Transportation forecasting package (originally INRO, Montreal; now Bentley OpenPaths). Strengths:
SOLA equilibrium assignment (faster-converging second-order linear approximation), extensive
transit assignment (frequency- and schedule-based), and a mature **Python API** (Modeller,
Desktop, Network, Matrix, Database APIs) that makes EMME among the most scriptable closed
packages; model-as-code workflows with JSON tool specifications. ODME implemented as a
generalised least-squares minimisation over route-choice probabilities (Bentley OpenPaths help).
Licence: no longer sold individually; bundled in OpenPaths Advanced (from ~$9,951/yr) or
Ultimate (from ~$12,877/yr, adds DYNAMEQ dynamic traffic assignment and AGENT).

### TransCAD (Caliper)
GIS-based travel demand modelling package (Caliper Corporation, USA), described as the most
widely used forecasting tool in US metropolitan planning (70%+ of MPOs per Caliper). Strengths:
tightly integrated GIS engine, fast multi-threaded UE assignment, dedicated **screenline
analysis** and assignment-difference utilities, flexible OD matrix estimation from counts, and
full demand suite (4-step, freight, activity models). Scripting via GISDK (Caliper Script).
Pricing is published: US$15,000 standard single workstation or US$6,000–7,200 annual; academic
licences restricted to unfunded research/teaching.

### CUBE (Bentley OpenPaths CUBE)
Land-use and transportation modelling suite (formerly Citilabs), organised around the Voyager
scripting language (NETWORK, MATRIX, HIGHWAY, PT, GENERATION, DISTRIBUTION, FRATAR programs)
with a flow-chart Application Manager, GIS integration via ESRI technology, CUBE Land for land
use, CUBE Avenue for mesoscopic, CubePy as a Python wrapper. Suited to large regional models
with congestion feedback loops between distribution and assignment. Licensing now via Bentley
OpenPaths bundles as above.

### OmniTRANS (Goudappel / Dat.mobility)
Dutch-standard multimodal modelling platform by Dat.mobility (part of Goudappel), used in ~80%
of Dutch municipal/regional traffic models per Goudappel. Modular: Essentials (networks, static
allocation, analysis), Matrix (OD estimation), Public Transport, Flow (STAQ semi-dynamic and
STREAMLINE dynamic car assignment with realistic congestion). Python integration and Octavius
person-level demand modelling. **Important transitional fact:** from 1 January 2026 OmniTRANS
Expert is superseded by "OmniTRANS Powered by Bentley OpenPaths" (OmniTRANS Expert maintained
until 2030), so the closed-tool market is consolidating further around Bentley's stack. Pricing
on request; free demo limited to 25 zones for 60 days.

---

## 3. What open tools are missing to fully replace SATURN-type workflows

1. **Junction-level simulation and its feedback on assignment.** SATURN's defining feature —
   cyclic flow profiles at signals, gap-acceptance-derived capacities (ARCADY/PICADY/OSCADY
   lineage), signal optimisation, and explicitly modelled blocking back/flow metering between
   adjacent junctions — has no open counterpart. AequilibraE and MATSim model congestion on
   links/queues but not at the level of individual junction control, lane structure and turning
   flows. This is the largest single capability gap for a Leeds-style urban model.
2. **Mature, validated ODME.** AequilibraE's matrix estimation exists but lacks the decades of
   practitioner track record behind SATME2, TFlowFuzzy or TransCAD's procedures; documented
   real-world UK applications are scarce.
3. **Bundled demand modelling.** Open stacks require assembling trip generation, distribution
   (partially available via AequilibraE's gravity/IPF) and mode choice yourself; closed packages
   ship integrated, convergent variable-demand loops.
4. **Certified validation tooling.** GEH calculations, screenline/cordon reporting, convergence
   audit trails and appraisal-grade output formats (e.g. integration with UK TAG appraisal
   software) are out-of-the-box in closed tools; in open tools they must be scripted.
5. **Institutional acceptance.** UK scheme appraisal, statutory consultation and peer review are
   calibrated around SATURN/standard commercial packages; open-model results face an
   evidential, not just technical, barrier.
6. **Support and continuity.** Commercial licences carry support contracts, training and
   indemnities; open tools rely on community/consultancy support.
7. **Parameter-ready networks from OSM.** Capacities, saturation flows, signal timings and HGV
   percentages are not in OSM; any OSM-based motor model still needs a substantial
   count/attribute fusion layer (osmactive-style OSM classification only covers active travel
   infrastructure).

## 4. Unique strengths of the open stack (and of SATURN)

**SATURN uniquely offers:** the cyclic-flow-profile junction simulation coupled to assignment;
blocking back/flow metering; UK-standard status with DfT/TfL/National Highways lineage; SATME2
with a long track record; Origin-Based Assignment for high-convergence multi-user-class runs;
SatGPU acceleration. No open tool currently reproduces this combination.

**Open tools uniquely offer:** zero licence cost and no licence servers; full source-code
accessibility and inspectability of algorithms; native headless/CI execution (AequilibraE CLI,
MATSim, sDNA CLI) enabling fully reproducible, version-controlled model pipelines; frictionless
integration with the Python/R data ecosystem (GeoPandas, Jupyter, OMX/SQLite data formats);
OSM-native network creation with parametrisable tag handling (AequilibraE, osmactive);
specialised strengths outside the SATURN frame — MATSim's fully dynamic agent-based behaviour
in one open framework, and cityseer/sDNA/pandana's pedestrian-scale centrality/accessibility
and flow heuristics that closed strategic packages do not cover. AequilibraE itself is the
closest open analogue to the SATURN assignment core: multi-class UE (FW/CFW/BFW/MSA) with
BPR-family VDFs, PCE/VoT user classes (e.g. a truck class with PCE 2.5), skimming, select-link
analysis and matrix estimation in an MIT-style scriptable pipeline.

---

## 5. Data-quality caveats

- **netassign / netweld / osmactive (AT Futures):** the task brief lists these as AT Futures
  tools. Web and GitHub searches could not locate tools named *netassign* or *netweld* under the
  ATFutures GitHub organisation (repos geoplumber, calendar, sumor, who, geobench, nyped found;
  the two candidate names returned 404s) or anywhere else on the public web. **osmactive** does
  exist, but as an R package from the NPT Scotland project (`nptscot/osmactive`) — not published
  under ATFutures — providing OSM extraction of drive/cycling networks and cycling level-of-service
  classification; it is not a motor assignment engine. The comparison above therefore covers
  osmactive as found and flags the other two names as unverifiable rather than describing them.
- SATURN's public help pages describe release 11.5 (2020) as "latest"; a 11.6 series and current
  2025 user-group activity are evidenced elsewhere, so feature details are dated to that
  documentation generation.
- Vendor pricing (PTV, Bentley, Caliper) is quoted from vendor/CIVITAS pages and changes; treat
  as indicative only.
- AequilibraE's licence: the repository metadata lists a non-standard licence identifier; confirm
  the exact licence text before redistribution.

## Sources

- SATURN help, "The Function of SATURN": https://saturnsoftware2.co.uk/SATURNhelp/Modules/31/Documents/thefunctionofsaturn.htm
- SATURN help, "Distribution of SATURN" (licensing): https://saturnsoftware2.co.uk/SATURNhelp/Modules/31/Documents/distributionofsaturn.htm
- SATURN product site: https://saturnsoftware2.co.uk/
- SATURN archive / User Group Meeting presentations: https://saturnsoftware2.co.uk/archive
- Saturn (software), Wikipedia: https://en.wikipedia.org/wiki/Saturn_(software)
- Using SATURN for CBD network modelling (GEH, SATME2, cyclic flow profiles): http://hdl.handle.net/10539/15739
- PTV Visum FAQs (assignment methods, ODME, licences): https://www.ptvgroup.com/en/products/ptv-visum/faqs
- PTV Visum product page: https://www.ptvgroup.com/en-us/products/ptv-visum
- CIVITAS on PTV Visum (price from ~€9,000): http://civitas.eu/solutions-and-tools/ptv-visum
- PTV Visum application examples (OSM import, TFlowFuzzy, SBA, DUE): https://cgi.ptvgroup.com/vision-help/VISUM_2025_ENG/Content/2_Herzlich_Willkommen/Examples_Anwendungsbeispiele.htm
- Bentley OpenPaths product/pricing: https://www.bentley.com/en/products/openpaths/
- Bentley OpenPaths EMME documentation: https://docs.bentley.com/LiveContent/web/OpenPaths-v2026/Help/en/topics/EMME/What_is_OpenPaths_EMME.html
- Bentley OpenPaths EMME Modeller API guide: https://docs.bentley.com/LiveContent/web/OpenPaths-v2026/Help/en/topics/EMME/Modeller/scripting_modeller.html
- Bentley OpenPaths matrix estimation (ODME maths): https://docs.bentley.com/LiveContent/web/OpenPaths-v2025/Help/en/topics/2653286/GUID-E7604B4F-98EF-4C10-9F38-90B98E93E471.html
- Bentley KB, "8 Steps to Building a CUBE Model" (Voyager programs): https://bentleysystems.service-now.com/community?id=kb_article&sysparm_article=KB0102094
- Bentley OpenPaths CubePy: https://docs.bentley.com/LiveContent/web/OpenPaths-v2025/Help/en/topics/2681859/GUID-0908905A-3943-4251-A3C0-26AD80E061E7.html
- Caliper TransCAD product page (assignment, screenline analysis, OD estimation): https://www.caliper.com/transcad/default.htm
- Caliper TransCAD pricing: https://www.caliper.com/tcprice.htm
- Caliper TransCAD travel demand brochure: https://www.caliper.com/pdfs/TravelDemandModelingBrochure.pdf
- Goudappel OmniTRANS Expert (modules, STAQ/STREAMLINE): https://www.goudappel.com/expertise/data-and-it-solutions/traffic-modelling-software-OmniTRANS-Expert
- Goudappel OmniTRANS Powered by Bentley OpenPaths: https://www.goudappel.com/expertise/Data-and-IT-solutions/traffic-modellingsoftware-OmniTRANS-Powered-by-Bentley-OpenPaths
- AequilibraE docs, Traffic Assignment Procedure: https://www.aequilibrae.com/develop/python/traffic_assignment/assignment_procedures.html
- AequilibraE docs, Traffic Assignment Insights (algorithms, convergence): https://www.aequilibrae.com/latest/python/traffic_assignment/traffic_assignment_insights.html
- AequilibraE docs, network import/export (OSM, GMNS): https://www.aequilibrae.com/latest/python/network_manipulation/import_and_export_network.html
- AequilibraE docs, Parameters YAML (OSM field mapping): https://www.aequilibrae.com/latest/python/aequilibrae_project/parameter_file.html
- AequilibraE docs, route choice & select link analysis: https://www.aequilibrae.com/docs/python/V.1.1.3/_auto_examples/assignment_workflows/plot_route_choice_basics.html
- QAequilibraE (QGIS GUI) traffic assignment and trip distribution menus: https://www.aequilibrae.com/latest/qgis/menus_in_detail/traffic_assignment.html and https://www.aequilibrae.com/develop/qgis/menus_in_detail/trip_distribution.html
- AequilibraE GitHub (releases incl. v1.7.0 CLI, licence note): https://github.com/AequilibraE/aequilibrae and https://github.com/AequilibraE/aequilibrae/releases
- MATSim: https://matsim.org/ and The MATSim book (part I): https://matsim.org/files/book/partOne-latest.pdf
- cityseer (GitHub / PyPI): https://github.com/benchmark-urbanism/cityseer-api and https://pypi.org/project/cityseer/
- sDNA open-source manual: https://sdna-open.readthedocs.io/en/latest/ ; Cardiff sDNA software page: https://research-sites.cardiff.ac.uk/spatial-design-network-analysis-sdna/software/
- sDNA+ GitHub (AGPL-3.0): https://github.com/fiftysevendegreesofrad/sdna_plus/
- pandana docs: http://udst.github.io/pandana/
- osmactive (NPT Scotland): https://github.com/nptscot/osmactive and https://nptscot.github.io/osmactive/
- ATFutures GitHub organisation (repos searched for netassign/netweld): https://github.com/ATFutures
- Practitioner context on 4-step modelling and GEH targets (arterials.co guide): https://www.arterials.co/how-to-build-a-transportation-travel-demand-model-a-practitioners-guide/

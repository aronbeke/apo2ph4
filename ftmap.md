# Apo2ph4 from FTMap/ATLAS hotspots

This document describes `apo2ph4_ftmap.py`, a driver that generates *apo* site
pharmacophores from [FTMap](https://ftmap.bu.edu/)/ATLAS output instead of from
docked fragments, and explains the design decisions, the caveats, and how to use
the resulting models downstream.

Nothing in the upstream Apo2ph4 workflow is modified. Two files are added:

| File | Role |
|---|---|
| `apo2ph4_ftmap.py` | driver: hotspot selection → grids → feature counts → interaction pharmacophores → final model |
| `scripts/ftmap_gen_ia_ph4s.py` | probe pharmacophore generation (default engine; also an in-memory replacement for `scripts/gen_ia_ph4s.py`) |
| `scripts/ph4_to_pymol.py` | visualise a `.pml` model in PyMOL (§5.4) |

---

## 1. Rationale

The stock Apo2ph4 workflow is four steps:

1. `apo2ph4_define_binding_site.py` — place a dummy ligand that defines the pocket.
2. `apo2ph4_prepare_and_dock.sh` — AutoGrid maps + Vina docking of a fragment
   library + pseudo co-crystal complexes + `feature_count.txt`.
3. `apo2ph4_generate_docked_frag_ph4s_cdpkit.sh` — one interaction pharmacophore
   per docked pose, concatenated into `pharmacophores/pharmacophores.pml`.
4. `apo2ph4_generate_ph4.py` — grid-guided scoring, clustering and selection of
   the features across all poses → the final model.

Step 4 is the method. Its inputs are (a) a multi-model PML of per-pose
interaction pharmacophores and (b) AutoGrid maps of the receptor. Step 2 exists
only to produce posed fragments in the pocket plus those maps.

FTMap already provides posed probe molecules, grouped into consensus clusters
that mark experimentally meaningful hotspots. Substituting them for the docked
fragments keeps steps 3 and 4 intact and removes the docking entirely:

```
FTMap consensus cluster SDF ─┐
                             ├─► pseudo-complexes ─► interaction ph4s ─┐
receptor PDB ────────────────┘                                         ├─► apo2ph4_generate_ph4.py ─► model.pml
                             └─► AutoGrid maps + feature_count.txt ────┘
```

Conceptually this is a better fit than docking: FTMap probes are placed by
mapping solvent-sized molecules onto the protein and retaining consensus sites,
which is exactly the hotspot information the Apo2ph4 scoring function tries to
recover from docked fragment density.

---

## 2. What the driver does, per site

For each site the driver creates `<out-dir>/site_<NN>_cs<ranks>/` containing:

```
probes.sdf              probe molecules of the site
receptor.pdb            receptor copy used for grids and complexes
gridbox_ligand.pdb      benzene translated to the probe centroid (grid box definition)
receptor.pdbqt          prepare_receptor4.py
receptor.gpf            prepare_gpf4.py
receptor.*.map          autogrid4  (C, e, HD, OA, A, d, ...)
receptor.glg            AutoGrid log
feature_count.txt       scripts/count_features.py over probes.sdf
complexes/              pseudo co-crystal PDBs (pdb mode, removed unless --keep-complexes)
complex_list.list       input list for gen_ia_ph4s.py
pharmacophores/
    pharmacophores.pml  one interaction pharmacophore per probe
    site_<NN>_cs<...>.pml   THE FINAL MODEL
ph4_report.txt          stdout of apo2ph4_generate_ph4.py (feature table)
site_summary.json       machine readable summary
```

and `<out-dir>/sites.json` indexing all of them.

### 2.1 Hotspot selection

FTMap consensus clusters are **separate sites**, not parts of one pocket. For
the bundled 1H24 example the 12 clusters span 13–34 Å:

```
CS0  (22 probes)  17.19  82.52  70.23    strongest hotspot
CS1  (14)         23.30  65.50  73.89
CS2  (12)         27.79  66.65  72.68    CS1/CS2/CS4 are 5–8 Å apart -> one druggable site
CS3  (11)         19.88  79.26  57.46
CS4   (9)         26.24  68.95  76.07
...
CS11  (2)         16.61  58.08  71.51
```

Running them together would be wrong twice over: a single 30 Å AutoGrid box
cannot cover them, and the density scoring in `apo2ph4_generate_ph4.py` compares
features within a 3 Å radius, so features from unrelated sites would pollute the
per-type normalisation and the ranking.

Therefore the driver processes one site at a time:

* `--clusters top` (default) — the highest ranked cluster only.
* `--clusters top:3` — the three highest ranked clusters.
* `--clusters 1,2,4` — explicit ranks.
* `--clusters all` — every cluster, one model each.

`--group-within D` (default 8 Å) merges selected clusters whose centroids lie
within `D` of each other into a single site, by single linkage. This is how
neighbouring consensus clusters — the classic FTMap signature of a druggable
pocket — are modelled as one pharmacophore rather than three. Set `0` to disable.

### 2.2 Grids without docking

`apo2ph4_generate_ph4.py` reads five maps by fixed name: `receptor.C.map`
(hydrophobic), `receptor.e.map` (electrostatic, used for PI/NI),
`receptor.HD.map` (donor), `receptor.OA.map` (acceptor), `receptor.A.map`
(aromatic). The driver reproduces the grid half of
`apo2ph4_prepare_and_dock.sh` and stops there:

```
prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
prepare_ligand4.py   -l gridbox_ligand.pdb -o gridbox_ligand.pdbqt
prepare_gpf4.py -l gridbox_ligand.pdbqt -r receptor.pdbqt -y \
                -p npts='100,100,100' -p ligand_types='C, OA, HD, A' -p spacing='0.300'
autogrid4 -p receptor.gpf
```

The dummy ligand is `data/benzene.sdf` translated to the probe centroid of the
site, which is how `apo2ph4_define_binding_site.py` defines a box upstream;
`-y` then centres the grid on it. `--grid-center X Y Z` overrides this and sets
`gridcenter` explicitly. Defaults give 100 × 0.3 Å = **30 Å box**; the driver
warns if the probe extent of a site exceeds it (raise `--npts`, or `--spacing`
at the cost of resolution).

No Vina run is performed, and `--skip-grids` reuses existing maps across reruns,
which is where almost all of the runtime sits.

### 2.3 feature_count.txt

`apo2ph4_generate_ph4.py` divides each feature's density score by the total
number of features of that type in the input molecule set
(`apo2ph4_generate_ph4.py:403`). Upstream this is counted over the whole
fragment library. The driver runs `scripts/count_features.py` over **the probe
molecules of that site**, so the normalisation matches the molecules actually
fed in. Using a mismatched count file silently skews the relative ranking of
feature types, which is the single easiest way to get a plausible-looking but
wrong model.

### 2.4 Pseudo-complexes — the subtle part

`scripts/gen_ia_ph4s.py:46` (`splitComplex`) identifies the ligand as **the
first atom in the file plus every atom sharing its residue code and sequence
number**. Upstream this works only incidentally: `pdbqt_to_cocrystal_complex.py`
relies on PyMOL writing objects in load order, with the receptor loaded last.

The driver therefore writes each complex explicitly as:

1. probe atoms as `HETATM`, residue `LIG Z 999`, serials `1..N`,
2. `CONECT` records,
3. the receptor `ATOM` records, renumbered to start at `N+1`.

Two details matter:

* **Serial renumbering.** The receptor's own serials start at 1, so without
  renumbering the probe's `CONECT` records would point at protein atoms.
* **Bond orders through `CONECT`.** `gen_ia_ph4s.py` never perceives bond
  orders; it only calls `perceiveSSSR`, `setRingFlags`,
  `calcImplicitHydrogenCounts`, `perceiveHybridizationStates` and
  `setAromaticityFlags`. The orders have to come from the file. The driver
  encodes a bond of order *n* by listing the partner *n* times, which CDPKit
  reads back via `PDB_DEDUCE_BOND_ORDERS_FROM_CONECT_RECORDS`. Without this,
  benzene would be read as cyclohexane: no aromatic feature, and two implicit
  hydrogens added per carbon instead of one.

The FTMap probes carry only polar hydrogens; `makeHydrogenComplete` plus
`calcHydrogen3DCoordinates` in `gen_ia_ph4s.py` fills in the rest, so no
protonation step is needed — provided the bond orders survived, hence the point
above.

### 2.5 Pharmacophore generation modes — and why `ligand` is the default

This is the part where FTMap input genuinely differs from docked fragments, and
it was only apparent once the pipeline was run.

`--ph4-mode pdb` and `--ph4-mode interaction` both build **interaction**
pharmacophores: a probe feature is kept only if it has a geometrically valid
partner in the receptor. `pdb` does it through the pseudo-complex PDBs and the
unmodified upstream `scripts/gen_ia_ph4s.py`; `interaction` does the same in
memory via `scripts/ftmap_gen_ia_ph4s.py`, which reads the receptor once and
combines it with each probe:

```python
compl_mol.assign(probe)        # probe first
tag_as_ligand(compl_mol)       # LIG / 999 / chain Z
compl_mol += receptor
```

splitting by that residue identity rather than by "first atom in the file". The
two give the same result; `interaction` is faster and keeps SD-file bond orders
without a PDB round trip.

**Measured on CS0 of the bundled 1H24 example, both discard most of the signal:
12 of 47 probe features survive** (10 hydrophobic, 1 aromatic, 1 HBD, **0 HBA**),
which is too few for step 4 to work at all. Relaxing the criteria barely helps:

| H···acceptor max | kept features |
|---|---|
| 2.8 Å (CDPKit default) | 12 — AR 1, H 10, HBD 1 |
| 3.2 Å | 14 — AR 1, H 10, HBD 3 |
| 3.5 Å | 16 — AR 1, H 10, HBD 5 |
| 4.0 Å | 18 — AR 1, H 10, HBD 7 |
| 3.5 Å + acceptor angle 180° | 17 — AR 1, H 10, HBD 8, **HBA 1** |

The reason is structural, not a bug: FTMap places probes by mapping and
clustering an energy function, not by optimising H-bond geometry against the
receptor, so the poses sit slightly too far out and mis-angled for the strict
criteria (CDPKit default: H···A 1.2–2.8 Å, A–H–D ≥ 130°, acceptor angle ≤ 85°).
In CS0 the probe environment is largely hydrophobic — ILE213, LEU214, TRP217,
GLN254 — and the one near-donor, GLN254 NE2, sits 3.63 Å from the benzaldehyde
carbonyl oxygen, a long H-bond by any standard.

`--ph4-mode ligand` (**default**) therefore keeps every probe feature and takes
the exclusion volumes from the receptor environment. The probe positions *are*
the hotspot information, and step 4 re-filters everything by AutoGrid energies
anyway.

That alone is not enough, though, because of a second interaction with step 4's
algorithm. `process_vector_features()`
(`apo2ph4_generate_ph4.py:435`) clusters HBD/HBA features by their **origin** —
`position ± length × orientation` — using `DBSCAN(eps=0.01, min_samples=3)`.
That demands three features whose origins coincide to within 0.01 Å. Interaction
pharmacophores satisfy it because their orientation points at the partner, so
origins land on receptor atoms shared between fragments. Probe features oriented
along their own bonds never coincide, `np.max(labels)+1` is 0, and **the entire
HBD and HBA lists are silently returned empty.**

So in `ligand` mode each H-bond feature is anchored on the nearest complementary
receptor feature within `--hb-pair-radius` (default 4.5 Å): orientation points
at it, length is set to that distance, geometry to `VECTOR`. The origin that
step 4 computes then lands exactly on the receptor feature, restoring the
invariant its clustering depends on. This is the strict H-bond criterion
replaced by a pure distance criterion — the pairing that FTMap poses can
actually support. Features with no complementary receptor feature in range are
dropped and counted (4 of 47 for CS0).

One further detail: probe acceptors are generated with `SPHERE` geometry even
though they carry an orientation, and the PML writer then emits a `<point>`
element, losing the orientation. Step 4 calls `Pharm.getOrientation()`
unconditionally for HBA (`apo2ph4_generate_ph4.py:359`) and raises
`ItemNotFound`. Forcing `VECTOR` geometry, as the anchoring step does, fixes
this too.

Use `pdb` or `interaction` if you want strict interaction semantics and accept a
sparse model, optionally with `--hb-max-length 3.5`. Use `ligand` for a usable
model.

---

## 3. Usage

```bash
conda activate apo2ph4
```

> **Note:** `python` inside the `apo2ph4` environment is Python **2.7.3**
> (shipped by MGLTools). Always invoke `python3`. The driver uses
> `sys.executable` for the sub-scripts, so it propagates the correct interpreter
> automatically.

Strongest hotspot of the bundled example:

```bash
python3 apo2ph4_ftmap.py atlas_output/1h24
```

The CS1/CS2/CS4 ensemble as one site, model also written as CDPKit `.cdf`:

```bash
python3 apo2ph4_ftmap.py atlas_output/1h24 --clusters 1,2,4 --group-within 10 --formats pml,cdf
```

Every cluster separately, 8 features each, pseudo-complexes retained:

```bash
python3 apo2ph4_ftmap.py atlas_output/1h24 --clusters all --group-within 0 -n 8 --keep-complexes
```

See the planned steps without executing anything:

```bash
python3 apo2ph4_ftmap.py atlas_output/1h24 --clusters all --dry-run
```

All relevant `apo2ph4_generate_ph4.py` options are forwarded (`-n`, `-t`, `-H`,
`-D`, `-P`, `-N`, `--H_energy`, `--PI_energy`, `--NI_energy`), and `--ph4-args`
passes anything else verbatim, e.g.
`--ph4-args "--HBA_ONLY_SPHERES 1 --REMOVE_UNREALISIC_AR_ANGLES 0"`.

### 3.1 Result on the bundled 1H24 example

`python3 apo2ph4_ftmap.py atlas_output/1h24` on CS0 (22 probes, 47 probe
features, 4 dropped as unpairable) produces a 7-feature model in about a minute,
almost all of it AutoGrid:

```
Type	  X	  Y	  Z	Score	Buriedness(PSP)		Grid energy [kcal/mol]
H   	18.15	83.23	69.49	62.37	 3 			-0.539
AR  	18.15	83.23	69.49	39.22	 3 			-0.476
HBA 	17.43	81.28	69.16	38.03	 4 			-0.735
HBD 	17.47	81.16	69.04	15.87	 4 			-0.06
HBD 	15.00	80.44	71.53	15.30	 5 			-0.055
HBA 	14.86	80.65	71.26	13.11	 5 			-0.773
H   	15.44	81.68	71.50	11.58	 5 			-0.667
```

plus 10 exclusion volumes. Every feature has a negative grid energy and a
buriedness of 3–5, and the H-bond vectors point at real receptor features — the
second HBA's origin is `11.159 81.564 69.817`, which is GLN254 NE2.

Two things to note when reading such a model:

* **The hydrophobic and aromatic features share coordinates** (18.15/83.23/69.49).
  Both come from the benzene/phenol ring centre; upstream Apo2ph4 does the same.
* **HBD and HBA features appear in near-coincident pairs** — 0.17 Å and 0.37 Å
  apart here. Tracing them back to the probe atoms they came from:

  | feature | sits exactly on | other probe atoms within 0.35 Å |
  |---|---|---|
  | HBA 17.43/81.28/69.16 | *tert*-butanol O | acetonitrile O, acetamide O, isopropanol O, acetaldehyde O, phenol O |
  | HBD 17.47/81.16/69.04 | isopropanol O | acetamide O, acetone O, *tert*-butanol O, phenol O |
  | HBD 15.00/80.44/71.53 | ethanol O | methylamine N |
  | HBA 14.86/80.65/71.26 | methylamine N | ethanol O |

  Two effects stack here. Each hydroxyl (and the amine) is a donor *and* an
  acceptor on the same atom, so one probe alone already produces both features
  at identical coordinates. On top of that, six chemically different probes put
  an oxygen within a third of an Ångström of the same point — which is not noise
  but the hotspot signal itself: that position is where FTMap says a polar atom
  belongs.

  For screening, requiring both simultaneously is only satisfiable by a
  hydroxyl- or amine-like group. If that is what you want to find, keep both; if
  a pure acceptor would do, mark the donor optional or drop it.

### 3.2 All 12 consensus clusters, and what each site costs

`python3 apo2ph4_ftmap.py atlas_output/1h24 --clusters all --formats pml,cdf`

The 12 consensus clusters collapse into **8 sites** under the default 8 Å
centroid merge (CS1+CS2+CS4 become one site, CS7+CS10+CS11 another). Three
sites yield a model; five do not — see below.

| Site | Clusters | Probes | Probe features | Model | Wall time |
|---|---|---|---|---|---|
| `site_00_cs000` | CS0 | 22 | 47 (43 kept) | **7** — H, AR, HBA, HBD, HBD, HBA, H | 31.1 s |
| `site_01_cs001-002-004` | CS1, CS2, CS4 | 35 | 67 (62 kept) | **9** — H, HBA, HBA, H, H, HBD, HBA, HBA, HBD | 34.8 s |
| `site_02_cs003` | CS3 | 11 | 21 (18 kept) | **4** — AR, H, HBA, HBD | 29.0 s |
| `site_03_cs005` | CS5 | 9 | 16 | *failed* | 30.3 s |
| `site_04_cs006` | CS6 | 4 | 9 | *failed* | 28.2 s |
| `site_05_cs007-010-011` | CS7, CS10, CS11 | 7 | 11 (10 kept) | *failed* | 28.8 s |
| `site_06_cs008` | CS8 | 3 | 11 | *failed* | 27.9 s |
| `site_07_cs009` | CS9 | 3 | 6 | *failed* | 27.9 s |

Total 238 s for all 8 sites on 6 cores.

**Per-site cost is essentially constant at ~28–35 s, and it is almost entirely
AutoGrid.** Breaking a site down (`site_00`, the largest):

| Step | Time | Scales with |
|---|---|---|
| receptor prep + AutoGrid (5 maps, 100³ points) | 25.1 s | receptor size and `--npts`³ — *not* the number of probes |
| `count_features.py` | 0.3 s | probes |
| probe pharmacophore generation (`ligand` mode) | 4.0 s | probes (0.18 s each) |
| `apo2ph4_generate_ph4.py` | 1.7 s | probe features |

So the number of probes barely matters: `site_06` with 3 probes took 27.9 s and
`site_01` with 35 took 34.8 s. If you need the run faster, the lever is
`--npts`/`--spacing` (grid cost goes as the cube), not the hotspot selection.
Dropping to `--npts 60` cuts the grid step roughly 5-fold. Conversely,
`--skip-grids` on a re-run with changed pharmacophore settings reduces a site to
about 6 s.

Disk is the other constant: 40 MB of maps per site, 310 MB for all eight.

**Why five sites produce nothing.** All five fail the same way — an
`IndexError` at `apo2ph4_generate_ph4.py:619`, which averages the top four
features' scores without checking that four exist. Fewer than four features
survive scoring and selection, so the run dies. Note this is *not* simply a
matter of too few probe features going in: `site_03` fed in 16 features (H 5,
Ar 2, HBD 3, HBA 6) and still ended under four, because the origin-DBSCAN in
`process_vector_features` needs three coincident features to form a cluster, and
with 9 probes the same interaction is rarely sampled three times.

The driver catches this per site, reports it and carries on with the remaining
sites, rather than aborting the whole run. The practical consequence is a floor
on usable hotspots: **roughly 10 probes per site**, matching the rule that the
top few consensus clusters are the ones worth modelling anyway.

---

## 4. Expected differences from a docked-fragment run

These are method-level consequences, worth knowing before interpreting output.

**Far fewer molecules.** A normal run docks ~200 fragments, often with two poses
each. CS0 of the example has 22 probes yielding 47 features. The density score
at `apo2ph4_generate_ph4.py:400` sums a sigmoidal distance term over
neighbouring features, and the automatic feature-count cutoff keeps features
scoring above 20 % of the mean of the top four (`:628`). With an order of
magnitude fewer input features the automatic mode returns small models — 7
features for CS0. That is a reasonable size, but consider setting `-n`
explicitly and widening `-t` to 4 if you want more.

**Step 4 crashes on very sparse input.** `score_top_four_avg` is computed by
indexing `pharmacophore_list[0..3]` unconditionally (`:616`), so if fewer than
four features survive filtering, the run dies with `IndexError: list index out
of range` rather than reporting an empty model. This is an upstream limitation;
seeing it means the input pharmacophores were too sparse — check the feature
counts that `-v` prints just above the traceback.

**Different chemistry.** The probe set is 16 small solvents. Their ligand
pharmacophore features, as CDPKit perceives them, are what the whole method has
to work with — for CS0 of the example:

| code | formula | features per probe |
|---|---|---|
| BNZ benzene | C6H6 | H, AR |
| CHX cyclohexane | C6H12 | H |
| EHN ethane | C2H6 | H x2 |
| IPH phenol | C6OH6 | H, AR, HBD, HBA |
| HBX benzaldehyde | C7OH6 | H, AR, HBA |
| CCN acetonitrile | C2NH3 | H, HBA |
| EOH ethanol | C2OH6 | HBD, HBA |
| IPA isopropanol | C3OH8 | HBD, HBA |
| TBU *tert*-butanol | C4OH10 | HBD, HBA |
| NME methylamine | CNH5 | **PI**, HBD, HBA |
| ACM acetamide | C2NOH5 | HBD, HBA |
| URE urea | CN2OH4 | HBD x2, HBA |
| ACE acetaldehyde | C2OH4 | HBA |
| ACN acetone | C3OH6 | HBA |
| DMF dimethylformamide | C3NOH7 | HBA |
| dimethyl ether | C2OH6 | HBA |

Acceptors dominate (19 of 47 features for CS0), donors follow (11), hydrophobes
are few and small (10), aromatics come only from benzene, phenol and
benzaldehyde (5). **Methylamine is the one ionizable probe**, and CDPKit does
assign it a positive ionizable feature (2 of 47 for CS0), so `-P` can fire from
the probe side; nothing in the set is negative ionizable, so `-N` can only come
from the receptor.

**Expect coincident HBD/HBA features.** Every hydroxyl and amine probe — EOH,
IPA, TBU, IPH, NME — generates a donor *and* an acceptor feature on the same
atom, because a hydroxyl oxygen both donates its O–H and accepts on its lone
pairs. Step 4 scores the two types in completely independent lists, with
separate caps, so nothing notices that a donor and an acceptor ended up at the
same point. See §3.1 for how this shows up in the final model.

**Tighter spatial spread.** Consensus clusters are compact (1.6–7.2 Å extent in
the example) compared with fragments docked across a whole pocket. Features are
more concentrated, so scores are higher in absolute terms but less
discriminating between neighbours; the 3 Å `-t` threshold does more merging.

**Probe energies are not used.** `ATLAS_ENERGY`,
`ATLAS_PROBE_CLUSTER_BOLTZMANN_ENERGY` and the cluster sizes are recorded in
`site_summary.json` but do not weight the pharmacophore. Weighting features by
their probe's cluster population would be a natural extension — the SD tags are
parsed and available.

---

## 5. Output and downstream use

### 5.1 Formats

The final model is a LigandScout `.pml`, written by
`Pharm.PMLFeatureContainerWriter`. `--formats pml,cdf` additionally writes
CDPKit's native binary `.cdf`, which is lossless and faster to load. Any
extension CDPKit's `Pharm.FeatureContainerWriter` recognises can be requested;
the conversion happens in-process after generation.

### 5.2 Machine readable summary

`site_summary.json` per site, plus `sites.json` indexing them:

```json
{
  "site": "site_00_cs000",
  "consensus_clusters": [
    {"rank": 0, "size": 22, "centroid": [17.19, 82.52, 70.23],
     "best_atlas_energy": -13.69, "source_file": "..."}
  ],
  "n_probes": 22,
  "site_centroid": [17.19, 82.52, 70.23],
  "pharmacophore_pml": "/.../site_00_cs000.pml",
  "features": [
    {"type": "HBA", "x": 16.9, "y": 82.1, "z": 70.4,
     "score": 41.2, "buriedness_psp": 5, "grid_energy": -0.71}
  ]
}
```

The `features` array is parsed from the table that `apo2ph4_generate_ph4.py`
only *prints*: per-feature score, PSP buriedness (0–7, how enclosed the point
is) and grid energy. Upstream these numbers are discarded after being echoed to
the terminal; they are exactly what you want for ranking features, filtering
sites, or deciding which hotspot to pursue, so the driver persists them.

### 5.3 Molecular fitting

The CDPKit Python bindings in the `apo2ph4` environment (version 1.3.0) provide
the full screening stack:

| Task | API |
|---|---|
| build a screening database | `Pharm.PSDScreeningDBCreator` |
| screen conformers against the model | `Pharm.ScreeningProcessor` + `Pharm.PSDScreeningDBAccessor` |
| align a single molecule | `Pharm.PharmacophoreAlignment`, `Pharm.SpatialFeatureAlignment` |
| score a fit | `Pharm.PharmacophoreFitScore`, `Pharm.PharmacophoreFitScreeningScore` |
| collect hits | `Pharm.FileScreeningHitCollector` |

> The pip `CDPKit` wheel installs **bindings only** — there are no `psdscreen` /
> `psdcreate` command line tools in this environment. Screening must go through
> the Python API, or through LigandScout with the `.pml`.

A minimal fitting sketch:

```python
import CDPL.Pharm as Pharm, CDPL.Base as Base

query = Pharm.BasicPharmacophore()
Pharm.PMLPharmacophoreReader(Base.FileIOStream('site_00_cs000.pml', 'r')).read(query)

proc = Pharm.ScreeningProcessor(Pharm.PSDScreeningDBAccessor('conformers.psd'))
proc.setHitCallback(lambda hit, score: print(hit.getHitMoleculeIndex(), score) or True)
proc.searchDB(query)
```

Exclusion volumes are included in the model, so the fit is shape-aware; drop
them with `--no-exclusion-volumes` in `ftmap_gen_ia_ph4s.py` if a purely feature
based query is wanted.

### 5.4 Visual inspection

The upstream README's advice applies unchanged, and more so here: **look at the
model together with the receptor** and confirm no feature sits outside the
pocket. With FTMap input there is a second thing to check — that the features
sit on the consensus cluster the site was built from — which is why the probes
are worth loading alongside.

`scripts/ph4_to_pymol.py` builds that scene in PyMOL, which is already in the
`apo2ph4` environment (LigandScout is not, and is the only thing that reads
`.pml` natively):

```bash
python3 scripts/ph4_to_pymol.py \
    ftmap_ph4/site_00_cs000/pharmacophores/site_00_cs000.pml \
    -r ftmap_ph4/site_00_cs000/receptor.pdb \
    -p ftmap_ph4/site_00_cs000/probes.sdf \
    -o ftmap_ph4/site_00_cs000/pharmacophores/site_00_cs000.pse
```

```bash
python3 -m pymol ftmap_ph4/site_00_cs000/pharmacophores/site_00_cs000.pse
```

> **Do not use the `pymol` command in this environment.** Its launcher script is
> `exec "$CONDA_PREFIX/bin/python" .../pymol/__init__.py`, and `python` here is
> the Python 2.7.3 that MGLTools pulls in, so it dies with a `SyntaxError` in
> PyMOL's own startup code. `python3 -m pymol` starts the same GUI correctly.
> On a headless shell, set `DISPLAY` first (this machine uses `:1`). Start it
> from a terminal with a tty: PyMOL reads commands from stdin, so launched as a
> detached background process it hits EOF and exits cleanly the moment the
> session finishes loading.

`--show` opens the GUI directly instead of writing a session, and `--png FILE`
renders an image headlessly (useful over SSH or in a batch loop over sites).

`--hide-probes` gives a clean view of the model alone: the probes are still
loaded, but switched off, so they stay in the object panel and one click brings
the cluster back. Omitting `-p` entirely leaves them out altogether.

What it draws:

| Element | Rendering |
|---|---|
| hydrophobic | yellow sphere |
| aromatic | blue sphere + a disc in the ring plane, both faces |
| H-bond acceptor | red sphere + arrow to the receptor donor it was paired with |
| H-bond donor | green sphere + arrow to the receptor acceptor |
| positive / negative ionizable | blue / magenta sphere |
| exclusion volume | translucent grey sphere (`--xvol-as-mesh` for a mesh) |
| receptor | translucent cartoon, side chains within `--env-radius` as sticks |
| probes | grey lines, all molecules at once |

The cartoon traces only the backbone, so the side-chain sticks read as a scatter
of small disconnected fragments around the ribbon. They are protein — the
receptor PDB contains nothing but standard amino acids — and they are what the
H-bond arrows point at, since the features were anchored on receptor donor and
acceptor atoms. With the defaults, 27 residues (117 atoms) are drawn this way
for CS0; lower `--env-radius` to thin them out.

Sphere radii are the feature tolerances, so what you see is the actual search
volume. Features are **pseudoatoms**, not flat CGO, so they stay clickable and
measurable — `distance` between a feature and a side-chain atom works as usual,
and each carries a residue name (`H`, `HBA`, …) and index for selection. The
arrows are the one exception; they are CGO and purely decorative.

Useful checks once the scene is open:

* Do the arrows end on chemically sensible receptor atoms? In `ligand` mode they
  were assigned by distance alone (§2.5), so this is where you catch a donor
  pointed at something it could never reach.
* Do coincident HBD/HBA pairs correspond to one hydroxyl, or to two genuinely
  different probes?
* Is any feature sticking out into solvent rather than into the pocket?

Other viewers: `.pml` is LigandScout's own format, so LigandScout gives the
richest view if you have a license. Chimera and VMD can show the receptor and
probes but have no notion of pharmacophore features, so they would need the same
kind of conversion.

---

## 6. Known limitations

* **One receptor conformation.** FTMap was run on a single structure; the model
  inherits whatever side-chain arrangement that structure has. Apo2ph4's
  docking-based route has the same limitation.
* **Receptor protonation.** `*_atlas_protein.pdb` carries the hydrogens FTMap
  added. They are used as-is; no independent protonation or pH assignment is
  performed.
* **Grid box vs. merged sites.** Merging distant clusters with a large
  `--group-within` can produce a site wider than the AutoGrid box. The driver
  warns but does not auto-resize, because enlarging the box at fixed `--npts`
  would silently coarsen the grid.
* **No probe weighting.** See §4.
* **`ligand` mode pairs by distance only.** A probe acceptor is anchored on the
  nearest receptor donor within 4.5 Å regardless of angle, so a feature can be
  pointed at a partner it could not actually hydrogen bond to. The vector's
  direction is therefore a clustering device and a hint, not an assertion about
  geometry. Tighten `--hb-pair-radius` if this matters to you.
* **Output size.** Each site keeps its AutoGrid maps, ~40 MB per site. Add
  `ftmap_ph4/` to `.gitignore`, or delete `receptor.*.map` once the models are
  generated.
* **`--skip-grids` does not validate** that existing maps were generated for the
  same centre and box as the current run. Delete the site directory when
  changing `--npts`, `--spacing` or `--grid-center`.
* **Automatic feature count is unreliable** with few probes. Use `-n`.
* **Sites with fewer than four surviving features crash step 4.**
  `apo2ph4_generate_ph4.py:619` sums `pharmacophore_list[0..3]` unconditionally.
  This is upstream behaviour, left untouched; the driver catches the failure,
  reports the site as skipped and continues. On the 1H24 example this costs
  five of eight sites (§3.2). A one-line upstream guard — dividing by
  `min(4, len(pharmacophore_list))` — would let these emit 1–3 feature models,
  though a three-feature pharmacophore is weak for screening anyway.

---

## 7. Reference: the example dataset

`atlas_output/1h24/` holds a standard ATLAS deliverable:

| File | Content |
|---|---|
| `*_consensus.<rank>.<size>.sdf` | 12 consensus clusters, 22→2 probes each, with `ATLAS_ENERGY`, `ATLAS_PROBE_CLUSTER_SIZE`, `ATLAS_CONSENSUS_CLUSTER_RANK` SD tags |
| `*_atlas_protein.pdb` | the protein FTMap used (hydrogens added, B-factors zeroed) — the receptor the driver picks by default |
| `*_atlas.pdb` | protein + all probes as HETATM residues + `LIG` atoms marking cluster centres |
| `*_chainR.pdb` | the FTMap input: chain R, ATOM records only, genuinely apo |
| `*_atlas.tar.xz` | the full FTMap archive |

The driver auto-detects `*_atlas_protein.pdb`; `--receptor` overrides it.

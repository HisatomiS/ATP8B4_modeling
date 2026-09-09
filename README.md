# ATP8B4-CDC50A Homology Model (MODELLER)

## Target
- Chain A: ATP8B4 (probable phospholipid-transporting ATPase IM), NP_079113.2, Homo sapiens
- Chain B: CDC50A (TMEM30A), NP_060717.1, Homo sapiens (identical sequence to template chain B)

## Template
- PDB: 8OXC (cryo-EM, resolution 2.58 Å)
- Chain A: ATP8B1 (AT8B1_HUMAN, PI-bound state)
- Chain B: CDC50A (CC50A_HUMAN)
- Ligands present in template (not carried into this model): BMA, MG, NAG, PIE (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoinositol), VN4 (vanadate)
- Template structure resolved from residue 63 of chain A (N-terminal 1-62 unresolved)

## Environment
- MODELLER 10.8, r13157
- WSL2 (Ubuntu), conda env `scanpy`, Python 3.11
- License key set via `export KEY_MODELLER='...'`

## Files
- `target.ali` — target sequence (ATP8B4 + CDC50A), PIR format, multi-chain
- `template.pdb` — 8OXC coordinates (chains A, B only; no HETATM)
- `01_get_template_seq.py` — extracts resolved sequence from template.pdb (`8oxc.seq`)
- `02_align.py` — align2d alignment of target vs. template (`alignment.ali`, `alignment.pap`)
- `03_model.py` — AutoModel run: 10 models, `md_level=refine.slow`, `repeat_optimization=2`
- `04_evaluate.py` — DOPE-based ranking of the 10 models
- `05_annotate_dope_profile.py` — maps template-uncovered (ab initio) residue ranges onto the DOPE profile

## Pipeline (run in order)
```bash
conda activate scanpy
export KEY_MODELLER='MODELIRANJE'

python 01_get_template_seq.py
python 02_align.py
python 03_model.py
python 04_evaluate.py
python 05_annotate_dope_profile.py
```

## Alignment summary
- TM domains and catalytic core: high identity, no significant gaps
- Regions with no template coverage (built ab initio, PI4-ATPase numbering, target/ungapped residue index):
  - N-terminal ~1-10 (ATP8B4 N-terminus)
  - ~420-430
  - ~700-760 (cytoplasmic loop near N-domain)
  - ~1040-1060
  - ~1100-1220 (C-terminal cytoplasmic tail)
  - ~1545-1553 (CDC50A N-terminal region)

## Model results
10 models generated (`ATP8B4_CDC50A.B99990001.pdb` – `B99990010.pdb`), ranked by DOPE score:

| Model | DOPE score |
|---|---|
| B99990001 | -185184.80 |
| B99990009 | -184897.31 |
| B99990010 | -184887.11 |
| B99990007 | -184610.73 |
| B99990008 | -184318.52 |
| B99990002 | -183976.27 |
| B99990005 | -183795.11 |
| B99990003 | -183064.91 |
| B99990006 | -182878.81 |
| B99990004 | -182480.58 |

**Best model: `ATP8B4_CDC50A.B99990001.pdb`** (DOPE = -185184.80)

DOPE profile inspection confirmed that local energy peaks correspond to the ab-initio (template-uncovered) regions listed above; no unexpected peaks were found in template-covered regions overlapping known mutation sites of interest.

## Variant: G395S

### Files
- `06_mutate_G395S.py` — standard MODELLER mutate_model.py workflow; mutates
  chain A residue 395 (G→S) on the best WT model and locally re-optimizes
  (MD + conjugate gradients) → `ATP8B4_CDC50A.B99990001SER395.pdb`
- `08_check_phi_psi.py` — backbone phi/psi angles around residue 395 (Biopython)
- `09_spatial_neighbors_G395.py` — residues within 8 Å of G395 in 3D,
  independent of sequence position (Biopython NeighborSearch)
- `10_check_g395s_clash.py` — phi angle before/after mutation; side-chain
  distance check (CB, OG) vs. neighbor residues
- `11_rotamer_scan_G395S.py` — full 360° chi1 scan of the Ser395 side chain
  (OG atom) vs. Asp815; separately checks whether the Cβ position (fixed
  by backbone, independent of chi1) is in contact
- `14_wt_vs_mut_asp815_distance.py` — WT Gly395(CA)–Asp815 distance;
  G395S Ser395 CB and OG distances to Asp815, reported separately

### Results

**DOPE (whole model):** WT = -185184.80, G395S = -185202.11

**Backbone conformation:** WT G395 phi = +92.8°. G395S phi = +81.9°.

**3D spatial neighbors of G395 (independent of sequence position, ≤8 Å,
chain A unless noted):** 650, 651, 653, 813, 814, 815(Asp), 816, 819, 832,
833, 834, 835, 849, 850, 851, 169, 170, 368.

**Distances to Asp815 (G395S model):**
- Cβ(Ser395) – OD2(Asp815) = 2.36 Å
- OG(Ser395) – OD2(Asp815) = 2.62 Å

**WT reference distance:** CA(Gly395) – OD2(Asp815) = 3.45 Å

**Chi1 rotamer scan (Ser395 OG vs. Asp815, 15° steps, 360°):**
- Farthest: chi1 = -120°, distance = 2.91 Å
- Closest: chi1 = 60°, distance = 1.65 Å
- Cβ position does not change with chi1 (fixed by backbone geometry).

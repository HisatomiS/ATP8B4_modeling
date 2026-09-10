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
- .py files in Github 


## Pipeline (run in order)
```bash

python 01_get_template_seq.py
python 02_align.py
python 03_model.py 
python 04_evaluate.py
python 05_annotate_dope_profile.py
python 06_mutate_G395S.py
python 08_check_phi_psi.py
python 09_spatial_neighbors_G395.py
python 10_check_g395s_clash.py
python 11_rotamer_scan_G395S.py
python 14_wt_vs_mut_asp815_distance.py
python 15_map_atp8b4_to_atp8b1.py
python 18_find_motif_resnums.py
python 16_mutate_G457S_ATP8B1.py
python 19_spatial_neighbors_G457_ATP8B1.py
python 20_wt_vs_mut_asp893_distance_ATP8B1.py
python 21_find_all_motifs.py
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

## Control: same mutation on the ATP8B1 experimental structure

To rule out the clash being a homology-modeling artifact, the equivalent
G→S mutation was made directly on the ATP8B1 cryo-EM structure
(`template.pdb`, 8OXC) rather than the ATP8B4 homology model.

Residue numbers were located directly from the template sequence via
motif search (DKTGT, GDGAND), not by offset arithmetic from the
alignment — offset-based mapping gave incorrect numbers (444/825) and was
discarded.

### Files
- `18_find_motif_resnums.py` — locates the actual ATP8B1 residue numbers
  for the DKTGT and GDGAND motifs directly from `template.pdb`
- `16_mutate_G457S_ATP8B1.py` — same mutate_model.py workflow as 06, run
  on `template.pdb` chain A residue 457 (G→S) → `template_G457S.pdb`
- `19_spatial_neighbors_G457_ATP8B1.py` — residues within 8 Å of Gly457 in
  the WT template, independent of sequence position
- `20_wt_vs_mut_asp893_distance_ATP8B1.py` — WT Gly457(CA)–Asp893 distance;
  G457S Ser457 CB and OG distances to Asp893, reported separately

### Results
- DKTGT motif: D454-K455-T456-G457-T458 (ATP8B1 numbering)
- GDGAND motif: G892-D893-G894-A895-N896-D897 (ATP8B1 numbering)
- Spatial neighbors of Gly457 (WT template, ≤8 Å, chain A): 458, 456, 711,
  460, 893(Asp), 455, 462, 910, 459, 454, 927, 911, 913, 892, 233, 929,
  928, 463, 430, 709, 234, 708, 897(Asp), 912, 894, 891, 732, 926

**Distances to Asp893:**
- WT: CA(Gly457) – OD2(Asp893) = 3.47 Å
- G457S: CB(Ser457) – OD2(Asp893) = 2.47 Å
- G457S: OG(Ser457) – OD2(Asp893) = 2.66 Å

For comparison, the ATP8B4 homology model gave: WT CA–Asp815 = 3.45 Å,
G395S CB–Asp815 = 2.36 Å, G395S OG–Asp815 = 2.62 Å.

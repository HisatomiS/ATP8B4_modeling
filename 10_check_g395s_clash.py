import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

WT_PDB = "ATP8B4_CDC50A.B99990001.pdb"
MUT_PDB = "ATP8B4_CDC50A.B99990001SER395.pdb"
CHAIN_ID = "A"
TARGET_RESIDUE = 395

# 09 positive residue
WATCH_RESIDUES = [650, 651, 653, 813, 814, 815, 816, 819,
                   832, 833, 834, 835, 849, 850, 851]

CLASH_THRESHOLD_A = 2.4 

def get_phi(pdb_file, chain_id, resnum):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_file)
    chain = structure[0][chain_id]
    residues = {res.id[1]: res for res in chain if res.id[0] == " "}
    res = residues[resnum]
    prev_res = residues.get(resnum - 1)
    if prev_res is None:
        return None
    c_prev = prev_res["C"].get_vector()
    n = res["N"].get_vector()
    ca = res["CA"].get_vector()
    c = res["C"].get_vector()
    return np.degrees(calc_dihedral(c_prev, n, ca, c))


# --- (A) phi ---
phi_wt = get_phi(WT_PDB, CHAIN_ID, TARGET_RESIDUE)
phi_mut = get_phi(MUT_PDB, CHAIN_ID, TARGET_RESIDUE)

print("=== (A) change phi ===")
print(f"WT  (Gly395) phi = {phi_wt:.1f}")
print(f"MUT (Ser395) phi = {phi_mut:.1f}")
diff = phi_mut - phi_wt
print(f"{diff:.1f}")
if abs(diff) > 20:
    print("→ backborn affected by mutation")
else:
    print("→ no shift of backnorn")

# --- (B) clash check ---
parser = PDBParser(QUIET=True)
mut_structure = parser.get_structure("mut", MUT_PDB)
mut_chain = mut_structure[0][CHAIN_ID]
mut_residues = {res.id[1]: res for res in mut_chain if res.id[0] == " "}

ser395 = mut_residues[TARGET_RESIDUE]
if ser395.get_resname() != "SER":
    print(f"\n warning: not 395S(infact: {ser395.get_resname()})")

sidechain_atoms = [a for a in ser395 if a.get_name() in ("CB", "OG")]

print("\n=== (B) Ser395側鎖と周辺残基の距離チェック ===")
print(f"{'residuenumber':>8} {'aa':>8} {'chainatom':>8} {'interactionatom':>8} {'distal(A)':>8}  判定")

any_clash = False
for resnum in WATCH_RESIDUES:
    if resnum not in mut_residues:
        continue
    other_res = mut_residues[resnum]
    min_dist = None
    min_pair = None
    for sc_atom in sidechain_atoms:
        for other_atom in other_res:
            if other_atom.element == "H":
                continue
            d = sc_atom - other_atom
            if min_dist is None or d < min_dist:
                min_dist = d
                min_pair = (sc_atom.get_name(), other_atom.get_name())
    if min_dist is None:
        continue
    judge = "!! CLASH" if min_dist < CLASH_THRESHOLD_A else "OK"
    if min_dist < CLASH_THRESHOLD_A:
        any_clash = True
    print(f"{resnum:>8} {other_res.get_resname():>8} {min_pair[0]:>8} "
          f"{min_pair[1]:>8} {min_dist:>8.2f}  {judge}")

print()
if any_clash:
    print("→ clash positive")
else:
    print("→ no clash from threshold A2.4")
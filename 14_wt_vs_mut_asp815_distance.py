from Bio.PDB import PDBParser

WT_PDB = "ATP8B4_CDC50A.B99990001.pdb"
MUT_PDB = "ATP8B4_CDC50A.B99990001SER395.pdb"
CHAIN_ID = "A"
TARGET_RESIDUE = 395
ASP_RESIDUE = 815

parser = PDBParser(QUIET=True)


def get_residue(pdb_file, resnum):
    structure = parser.get_structure("s", pdb_file)
    chain = structure[0][CHAIN_ID]
    for res in chain:
        if res.id[0] == " " and res.id[1] == resnum:
            return res
    return None


def min_dist_to_residue(atom_list, other_res):
    asp_atoms = [a for a in other_res if a.element != "H"]
    best = None
    best_pair = None
    for atom in atom_list:
        for other in asp_atoms:
            d = atom - other
            if best is None or d < best:
                best = d
                best_pair = (atom.get_name(), other.get_name())
    return best, best_pair


# --- (1) WT: Gly395 (CA) - Asp815 ---
wt_res = get_residue(WT_PDB, TARGET_RESIDUE)
wt_asp = get_residue(WT_PDB, ASP_RESIDUE)
print(f"WT 395: {wt_res.get_resname()}")

wt_ca_atoms = [a for a in wt_res if a.get_name() == "CA"]
dist, pair = min_dist_to_residue(wt_ca_atoms, wt_asp)
print(f"WT: CA(395) - Asp815:{pair[1]}  distance = {dist:.2f} A")

# --- (2) MUT: Ser395 CB / OG -Asp815 ---
mut_res = get_residue(MUT_PDB, TARGET_RESIDUE)
mut_asp = get_residue(MUT_PDB, ASP_RESIDUE)
print(f"\nG395S : {mut_res.get_resname()}")

for atom_name in ("CB", "OG"):
    atoms = [a for a in mut_res if a.get_name() == atom_name]
    if not atoms:
        print(f"{atom_name}: no atom")
        continue
    dist, pair = min_dist_to_residue(atoms, mut_asp)
    print(f"MUT: {atom_name}(395) - Asp815:{pair[1]}  distance = {dist:.2f} A")
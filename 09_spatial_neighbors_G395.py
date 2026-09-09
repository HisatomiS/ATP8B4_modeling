from Bio.PDB import PDBParser, NeighborSearch

WT_PDB = "ATP8B4_CDC50A.B99990001.pdb"
CHAIN_ID = "A"
TARGET_RESIDUE = 395
CUTOFF_A = 8.0  # Å

parser = PDBParser(QUIET=True)
structure = parser.get_structure("wt", WT_PDB)
model = structure[0]

all_atoms = [atom for chain in model for atom in chain.get_atoms()]
ns = NeighborSearch(all_atoms)

target_chain = model[CHAIN_ID]
target_residue = None
for res in target_chain:
    if res.id[0] == " " and res.id[1] == TARGET_RESIDUE:
        target_residue = res
        break

if target_residue is None:
    raise SystemExit(f"residue{TARGET_RESIDUE}notfound")

nearby_residues = {}  # (chain_id, resnum) -> nearest

for atom in target_residue:
    close_atoms = ns.search(atom.coord, CUTOFF_A, level="A")
    for close_atom in close_atoms:
        parent_res = close_atom.get_parent()
        parent_chain = parent_res.get_parent()
        if parent_res.id[0] != " ":
            continue  
        key = (parent_chain.id, parent_res.id[1])
        if key == (CHAIN_ID, TARGET_RESIDUE):
            continue  
        dist = atom - close_atom
        if key not in nearby_residues or dist < nearby_residues[key][0]:
            nearby_residues[key] = (dist, parent_res.get_resname())

sorted_neighbors = sorted(nearby_residues.items(), key=lambda x: x[1][0])

print(f"fromG395within{CUTOFF_A}Å:")
print(f"{'chain':>4} {'residuenumber':>8} {'aa':>8} {'nearest(Å)':>12}")
for (chain_id, resnum), (dist, resname) in sorted_neighbors:
    seq_gap = abs(resnum - TARGET_RESIDUE) if chain_id == CHAIN_ID else None
    far_in_seq = ""
    if seq_gap is not None and seq_gap > 15:
        far_in_seq = "  <- confomational proximal"
    print(f"{chain_id:>4} {resnum:>8} {resname:>8} {dist:>12.2f}{far_in_seq}")
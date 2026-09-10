from Bio.PDB import PDBParser, NeighborSearch

WT_PDB = "template.pdb"
CHAIN_ID = "A"
TARGET_RESIDUE = 457
CUTOFF_A = 8.0

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
    raise SystemExit(f"residue {TARGET_RESIDUE} not found")

nearby_residues = {}

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

print(f"Residues within {CUTOFF_A} A of Gly{TARGET_RESIDUE}, sorted by distance:")
print(f"{'chain':>6} {'resnum':>8} {'AA':>5} {'min dist (A)':>14}")
for (chain_id, resnum), (dist, resname) in sorted_neighbors:
    print(f"{chain_id:>6} {resnum:>8} {resname:>5} {dist:>14.2f}")

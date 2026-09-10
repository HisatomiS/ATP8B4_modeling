from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one

WT_PDB = "template.pdb"
CHAIN_ID = "A"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("wt", WT_PDB)
chain = structure[0][CHAIN_ID]

residue_list = []
for res in chain:
    if res.id[0] != " ":
        continue
    try:
        one = index_to_one(three_to_index(res.get_resname()))
    except KeyError:
        one = "X"
    residue_list.append((res.id[1], one))

seq = "".join(r[1] for r in residue_list)


def find_motif(motif):
    idx = seq.find(motif)
    if idx == -1:
        print(f"motif '{motif}' not found")
        return
    print(f"motif '{motif}' found at sequence index {idx}")
    for i in range(idx, idx + len(motif)):
        resnum, aa = residue_list[i]
        print(f"  resnum {resnum}: {aa}")


find_motif("SDKTGT")
print()
find_motif("GDGAND")

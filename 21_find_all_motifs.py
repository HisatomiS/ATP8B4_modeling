from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one

TEMPLATE_PDB = "template.pdb"
TEMPLATE_CHAIN = "A"
TARGET_ALI = "target.ali"
TARGET_CODE = "ATP8B4_CDC50A"

MOTIFS = ["SDKTGT", "GDGAND", "TGDN"]


def get_template_residues():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("wt", TEMPLATE_PDB)
    chain = structure[0][TEMPLATE_CHAIN]
    residue_list = []
    for res in chain:
        if res.id[0] != " ":
            continue
        try:
            one = index_to_one(three_to_index(res.get_resname()))
        except KeyError:
            one = "X"
        residue_list.append((res.id[1], one))
    return residue_list


def get_target_residues():
    seqs = {}
    code = None
    buf = []
    with open(TARGET_ALI) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        if line.startswith(">P1;"):
            if code is not None:
                seqs[code] = "".join(buf).replace("/", "")
            code = line[4:].strip()
            buf = []
            i += 2
            continue
        else:
            buf.append(line.strip())
        i += 1
    if code is not None:
        seqs[code] = "".join(buf).replace("/", "")
    seq = seqs[TARGET_CODE].rstrip("*")
    # target.ali has no gaps for ATP8B4 itself (chain A is residues 1..N)
    residue_list = [(i + 1, aa) for i, aa in enumerate(seq)]
    return residue_list


def find_motif(residue_list, motif, label):
    seq = "".join(r[1] for r in residue_list)
    idx = seq.find(motif)
    if idx == -1:
        print(f"[{label}] motif '{motif}' not found")
        return
    resnums = [residue_list[i][0] for i in range(idx, idx + len(motif))]
    print(f"[{label}] '{motif}' -> residues {resnums[0]}-{resnums[-1]}")
    for i in range(idx, idx + len(motif)):
        resnum, aa = residue_list[i]
        print(f"    {resnum}: {aa}")


template_residues = get_template_residues()
target_residues = get_target_residues()

for motif in MOTIFS:
    find_motif(template_residues, motif, "ATP8B1 (template)")
    find_motif(target_residues, motif, "ATP8B4 (target)")
    print()

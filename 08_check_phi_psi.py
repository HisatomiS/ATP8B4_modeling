import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

WT_PDB = "ATP8B4_CDC50A.B99990001.pdb"
CHAIN_ID = "A"
CENTER_RESIDUE = 395
WINDOW = 5  # howmanyresidues

parser = PDBParser(QUIET=True)
structure = parser.get_structure("wt", WT_PDB)
chain = structure[0][CHAIN_ID]

residues = {res.id[1]: res for res in chain if res.id[0] == " "}

print(f"{'residue number':>8} {'aa':>8} {'phi':>10} {'psi':>10}")

for resnum in range(CENTER_RESIDUE - WINDOW, CENTER_RESIDUE + WINDOW + 1):
    if resnum not in residues:
        continue
    res = residues[resnum]
    resname = res.get_resname()

    prev_res = residues.get(resnum - 1)
    next_res = residues.get(resnum + 1)

    phi = psi = None
    try:
        if prev_res is not None:
            c_prev = prev_res["C"].get_vector()
            n = res["N"].get_vector()
            ca = res["CA"].get_vector()
            c = res["C"].get_vector()
            phi = np.degrees(calc_dihedral(c_prev, n, ca, c))
        if next_res is not None:
            n = res["N"].get_vector()
            ca = res["CA"].get_vector()
            c = res["C"].get_vector()
            n_next = next_res["N"].get_vector()
            psi = np.degrees(calc_dihedral(n, ca, c, n_next))
    except KeyError:
        pass

    comment = ""
    if resnum == CENTER_RESIDUE:
        comment += "<- G395 "
    if phi is not None and phi > 0:
        comment += "[phi > 0: suitable only for Gly]"

    phi_str = f"{phi:10.1f}" if phi is not None else f"{'--':>10}"
    psi_str = f"{psi:10.1f}" if psi is not None else f"{'--':>10}"
    print(f"{resnum:>8} {resname:>8} {phi_str} {psi_str}  {comment}")
import numpy as np 
from Bio.PDB import PDBParser 
 
MUT_PDB = "ATP8B4_CDC50A.B99990001SER395.pdb" 
CHAIN_ID = "A" 
TARGET_RESIDUE = 395 
ASP_RESIDUE = 815 
 
CLASH_THRESHOLD_A = 2.4 
 
 
def rotate_point_around_axis(point, axis_origin, axis_dir, angle_deg): 
    """Rodrigues' rotation formula: rotate around the axis by angle_deg degrees""" 
    angle = np.radians(angle_deg) 
    axis = axis_dir / np.linalg.norm(axis_dir) 
    p = point - axis_origin 
    p_rot = (p * np.cos(angle) 
             + np.cross(axis, p) * np.sin(angle) 
             + axis * np.dot(axis, p) * (1 - np.cos(angle))) 
    return p_rot + axis_origin 
 
 
def dihedral(p0, p1, p2, p3):
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.degrees(np.arctan2(y, x))
 
 
parser = PDBParser(QUIET=True) 
structure = parser.get_structure("mut", MUT_PDB) 
chain = structure[0][CHAIN_ID] 
residues = {res.id[1]: res for res in chain if res.id[0] == " "} 
 
ser395 = residues[TARGET_RESIDUE] 
asp815 = residues[ASP_RESIDUE] 
 
n = ser395["N"].get_coord() 
ca = ser395["CA"].get_coord() 
cb = ser395["CB"].get_coord() 
og = ser395["OG"].get_coord() 
 
asp_atoms = [(a.get_name(), a.get_coord()) for a in asp815 if a.element != "H"] 
 
# --- (A) Check for clashes between the CB atom itself and Asp815 --- 
print("=== Distance between the CB atom (backbone-derived, independent of chi1) and Asp815 ===") 
min_cb_dist = None 
min_cb_atom = None 
for name, coord in asp_atoms: 
    d = np.linalg.norm(cb - coord) 
    if min_cb_dist is None or d < min_cb_dist: 
        min_cb_dist = d 
        min_cb_atom = name 
print(f"CB - Asp815:{min_cb_atom}  distance = {min_cb_dist:.2f} A") 
if min_cb_dist < CLASH_THRESHOLD_A: 
    print("→ Clash at the CB level. Regardless of side-chain orientation, amino acids " 
          "(all amino acids except Gly) with a CB may not be able to fit at this position.") 
else: 
    print("→ No clash at the CB level.") 
 
# --- (B) Systematically vary chi1 from 0-360 degrees and examine the distance between OG and Asp815 --- 
print("\n=== Systematically vary chi1 (Ser side-chain orientation) from 0-360 degrees and calculate the minimum distance between OG and Asp815 ===") 
 
current_chi1 = dihedral(n, ca, cb, og) 
print(f"(Current chi1 in the model = {current_chi1:.1f} degrees)") 
 
axis_dir = cb - ca 
best_dist = -1 
best_angle = None 
worst_dist = 1e9 
worst_angle = None 
 
for target_chi1 in range(-180, 181, 15): 
    delta = target_chi1 - current_chi1 
    og_rot = rotate_point_around_axis(og, ca, axis_dir, delta) 
 
    min_dist = None 
    for name, coord in asp_atoms: 
        d = np.linalg.norm(og_rot - coord) 
        if min_dist is None or d < min_dist: 
            min_dist = d 
 
    if min_dist > best_dist: 
        best_dist = min_dist 
        best_angle = target_chi1 
    if min_dist < worst_dist: 
        worst_dist = min_dist 
        worst_angle = target_chi1 
 
print(f"chi1 with the greatest distance between OG and Asp815 = {best_angle} degrees (distance {best_dist:.2f} A)") 
print(f"chi1 with the smallest distance between OG and Asp815 = {worst_angle} degrees (distance {worst_dist:.2f} A)") 
 
print("\n=== Conclusion ===") 
if min_cb_dist < CLASH_THRESHOLD_A: 
    print("The CB atom itself clashes with Asp815, and this cannot be resolved by changing") 
    print("chi1 (side-chain orientation). This means that the clash is not simply due to") 
    print("an unfavorable orientation of the Ser side chain, but suggests that placing") 
    print("a residue with a CB (all amino acids except Gly) at this position is difficult") 
    print("at the backbone level.") 
else: 
    print("The CB atom itself does not clash with Asp815. The OG clash may be resolved") 
    print(f"depending on the choice of side-chain orientation (chi1={best_angle} degrees, distance={best_dist:.2f} A),") 
    print("suggesting that the optimization in STEP 06 may have been trapped in a local minimum.") 

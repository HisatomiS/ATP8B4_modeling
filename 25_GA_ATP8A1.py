import os
import re
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from modeller import *
from modeller.optimizers import MolecularDynamics
from modeller.automodel import autosched

PDB_FILE = "ATP8A1.pdb"
CHAIN_ID = "A"
OUT_BASENAME = "ATP8A1_GA"
NEIGHBOR_CUTOFF = 8.0
NEW_RESTYPE = "ALA"
NEW_ATOM_NAMES = ("CB",)

DKTGT_PATTERN = "DKTGT"
GDGAND_PATTERN = re.compile("GDG.ND")


def get_residue_list(pdb_file, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_file)
    chain = structure[0][chain_id]
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


def find_literal_motif(residue_list, motif):
    seq = "".join(r[1] for r in residue_list)
    idx = seq.find(motif)
    if idx == -1:
        return None
    return [residue_list[i] for i in range(idx, idx + len(motif))]


def find_regex_motif(residue_list, pattern):
    seq = "".join(r[1] for r in residue_list)
    m = pattern.search(seq)
    if m is None:
        return None
    start, end = m.span()
    return [residue_list[i] for i in range(start, end)]


residues = get_residue_list(PDB_FILE, CHAIN_ID)

dktgt = find_literal_motif(residues, DKTGT_PATTERN)
if dktgt is None:
    raise SystemExit(f"DKTGT motif not found in {PDB_FILE} chain {CHAIN_ID}")

print(f"DKTGT motif: {dktgt[0][0]}-{dktgt[-1][0]}")
for resnum, aa in dktgt:
    print(f"  {resnum}: {aa}")

gdgand = find_regex_motif(residues, GDGAND_PATTERN)
if gdgand is None:
    print("GDGAND-like motif not found")
else:
    print(f"\nGDGAND-like motif: {gdgand[0][0]}-{gdgand[-1][0]}")
    for resnum, aa in gdgand:
        print(f"  {resnum}: {aa}")

target_resnum = dktgt[3][0]
print(f"\nMutation target: Gly{target_resnum} -> {NEW_RESTYPE}")


def optimize(atmsel, sched):
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    atmsel.energy()


def refine(atmsel):
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                            md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                 (200, 600, (1000.0,))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                        max_iterations=its, equilibrate=equil)
            init_vel = False
    atmsel.energy()


def make_restraints(mdl1, aln):
    rsr = mdl1.restraints
    rsr.clear()
    s = Selection(mdl1)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(s, restraint_type=typ + '_dihedral', spline_range=4.0,
                  spline_dx=0.3, spline_min_points=5, aln=aln,
                  spline_on_site=True)


log.verbose()

env = Environ(rand_seed=-49837)
env.io.hetatm = False
env.edat.dynamic_sphere = False
env.edat.dynamic_lennard = True
env.edat.contact_shell = 4.0
env.edat.update_dynamic = 0.39

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl1 = Model(env, file=PDB_FILE)
ali = Alignment(env)
ali.append_model(mdl1, atom_files=PDB_FILE, align_codes=OUT_BASENAME)

s = Selection(mdl1.chains[CHAIN_ID].residues[str(target_resnum)])
s.mutate(residue_type=NEW_RESTYPE)

ali.append_model(mdl1, align_codes=OUT_BASENAME)

mdl1.clear_topology()
mdl1.generate_topology(ali[-1])
mdl1.transfer_xyz(ali)

mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

mdl2 = Model(env, file=PDB_FILE)
mdl1.res_num_from(mdl2, ali)

tmp_file = f"{OUT_BASENAME}.tmp"
mdl1.write(file=tmp_file)
mdl1.read(file=tmp_file)

make_restraints(mdl1, ali)
mdl1.env.edat.nonbonded_sel_atoms = 1

sched = autosched.loop.make_for_model(mdl1)

s = Selection(mdl1.chains[CHAIN_ID].residues[str(target_resnum)])
mdl1.restraints.unpick_all()
mdl1.restraints.pick(s)

s.energy()
s.randomize_xyz(deviation=4.0)

mdl1.env.edat.nonbonded_sel_atoms = 2
optimize(s, sched)

mdl1.env.edat.nonbonded_sel_atoms = 1
optimize(s, sched)

s.energy()

refine(s)
mdl1.env.edat.nonbonded_sel_atoms = 2
optimize(s, sched)

mut_file = f"{OUT_BASENAME}.pdb"
mdl1.write(file=mut_file)
os.remove(tmp_file)
print(f"\nWrote {mut_file}")

parser = PDBParser(QUIET=True)
wt_structure = parser.get_structure("wt", PDB_FILE)
wt_model = wt_structure[0]

all_atoms = [atom for chain in wt_model for atom in chain.get_atoms()]
ns = NeighborSearch(all_atoms)

wt_chain = wt_model[CHAIN_ID]
wt_target_res = None
for res in wt_chain:
    if res.id[0] == " " and res.id[1] == target_resnum:
        wt_target_res = res
        break

nearby = {}
for atom in wt_target_res:
    close_atoms = ns.search(atom.coord, NEIGHBOR_CUTOFF, level="A")
    for close_atom in close_atoms:
        parent_res = close_atom.get_parent()
        parent_chain = parent_res.get_parent()
        if parent_res.id[0] != " ":
            continue
        key = (parent_chain.id, parent_res.id[1])
        if key == (CHAIN_ID, target_resnum):
            continue
        dist = atom - close_atom
        if key not in nearby or dist < nearby[key][0]:
            nearby[key] = (dist, parent_res.get_resname())

sorted_neighbors = sorted(nearby.items(), key=lambda x: x[1][0])

print(f"\nResidues within {NEIGHBOR_CUTOFF} A of Gly{target_resnum} (WT), sorted by distance:")
print(f"{'chain':>6} {'resnum':>8} {'AA':>5} {'min dist (A)':>14}")
for (chain_id, resnum), (dist, resname) in sorted_neighbors:
    print(f"{chain_id:>6} {resnum:>8} {resname:>5} {dist:>14.2f}")

acidic_neighbors = [(c, r) for (c, r), (d, name) in sorted_neighbors if name in ("ASP", "GLU")]

mut_structure = parser.get_structure("mut", mut_file)
mut_chain = mut_structure[0][CHAIN_ID]
mut_residues = {res.id[1]: res for res in mut_chain if res.id[0] == " "}
wt_residues = {res.id[1]: res for res in wt_chain if res.id[0] == " "}

mut_target_res = mut_residues[target_resnum]

print(f"\nDistance check for nearby acidic (Asp/Glu) residues:")
for chain_id, resnum in acidic_neighbors:
    if chain_id != CHAIN_ID:
        continue
    wt_other = wt_residues[resnum]
    mut_other = mut_residues[resnum]

    wt_ca = [a for a in wt_target_res if a.get_name() == "CA"]
    wt_other_atoms = [a for a in wt_other if a.element != "H"]
    wt_dist = min(a - b for a in wt_ca for b in wt_other_atoms)

    print(f"\n{wt_other.get_resname()}{resnum}:")
    print(f"  WT CA(Gly{target_resnum}) - {wt_other.get_resname()}{resnum}  min distance = {wt_dist:.2f} A")

    for atom_name in NEW_ATOM_NAMES:
        mut_atoms = [a for a in mut_target_res if a.get_name() == atom_name]
        if not mut_atoms:
            continue
        mut_other_atoms = [a for a in mut_other if a.element != "H"]
        mut_dist = min(a - b for a in mut_atoms for b in mut_other_atoms)
        flag = "  !! < 2.4 A" if mut_dist < 2.4 else ""
        print(f"  MUT {atom_name}(Ala{target_resnum}) - {mut_other.get_resname()}{resnum}  min distance = {mut_dist:.2f} A{flag}")

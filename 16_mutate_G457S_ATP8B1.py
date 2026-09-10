import os
from modeller import *
from modeller.optimizers import MolecularDynamics
from modeller.automodel import autosched

MODEL_FILE = "template.pdb"
OUT_BASENAME = "template_G457S"
CHAIN = "A"
RESPOS = "457"
RESTYP = "SER"


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

mdl1 = Model(env, file=MODEL_FILE)
ali = Alignment(env)
ali.append_model(mdl1, atom_files=MODEL_FILE, align_codes=OUT_BASENAME)

s = Selection(mdl1.chains[CHAIN].residues[RESPOS])
s.mutate(residue_type=RESTYP)

ali.append_model(mdl1, align_codes=OUT_BASENAME)

mdl1.clear_topology()
mdl1.generate_topology(ali[-1])
mdl1.transfer_xyz(ali)

mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

mdl2 = Model(env, file=MODEL_FILE)
mdl1.res_num_from(mdl2, ali)

tmp_file = f"{OUT_BASENAME}.tmp"
mdl1.write(file=tmp_file)
mdl1.read(file=tmp_file)

make_restraints(mdl1, ali)
mdl1.env.edat.nonbonded_sel_atoms = 1

sched = autosched.loop.make_for_model(mdl1)

s = Selection(mdl1.chains[CHAIN].residues[RESPOS])
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

out_file = f"{OUT_BASENAME}.pdb"
mdl1.write(file=out_file)

os.remove(tmp_file)

print(f"Wrote {out_file}")

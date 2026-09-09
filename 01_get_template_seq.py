from modeller import *

TEMPLATE_PDB = "template.pdb"
TEMPLATE_CODE = "8oxc"

env = Environ()
env.io.hetatm = False
env.io.water = False

mdl = Model(env)
mdl.read(file=TEMPLATE_PDB, model_segment=('FIRST:A', 'LAST:B'))

aln = Alignment(env)
aln.append_model(mdl, align_codes=TEMPLATE_CODE, atom_files=TEMPLATE_PDB)
aln.write(file=f'{TEMPLATE_CODE}.seq')

print(f"Wrote {TEMPLATE_CODE}.seq")
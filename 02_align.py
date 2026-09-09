from modeller import *

TEMPLATE_PDB = "template.pdb"
TEMPLATE_CODE = "8oxc"
TARGET_CODE = "ATP8B4_CDC50A"

env = Environ()
env.io.hetatm = False
env.io.water = False

aln = Alignment(env)

mdl = Model(env, file=TEMPLATE_PDB, model_segment=('FIRST:A', 'LAST:B'))
aln.append_model(mdl, align_codes=TEMPLATE_CODE, atom_files=TEMPLATE_PDB)

aln.append(file='target.ali', align_codes=TARGET_CODE)

aln.align2d(max_gap_length=50)

aln.write(file='alignment.ali', alignment_format='PIR')
aln.write(file='alignment.pap', alignment_format='PAP')

print("done")

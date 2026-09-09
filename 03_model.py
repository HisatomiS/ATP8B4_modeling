from modeller import *
from modeller.automodel import *

env = Environ()
env.io.hetatm = False
env.io.water = False

a = AutoModel(env,
              alnfile='alignment.ali',
              knowns='8oxc',
              sequence='ATP8B4_CDC50A',
              assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 10
a.md_level = refine.slow
a.repeat_optimization = 2
a.max_molpdf = 1e6

a.make()

print("done")

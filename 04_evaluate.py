from modeller import *
from modeller.scripts import complete_pdb
import glob

env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

results = []
for pdb in sorted(glob.glob('ATP8B4_CDC50A.B*.pdb')):
    mdl = complete_pdb(env, pdb)
    s = Selection(mdl)
    dope_score = s.assess_dope(
        output='ENERGY_PROFILE NO_REPORT',
        file=pdb.replace('.pdb', '.profile'),
        normalize_profile=True,
        smoothing_window=15,
    )
    results.append((pdb, dope_score))

results.sort(key=lambda x: x[1])

print("ranking:")
for pdb, score in results:
    print(f"  {pdb}: DOPE = {score:.2f}")

best_pdb, best_score = results[0]
print(f"\n bestmodel: {best_pdb} (DOPE = {best_score:.2f})")
print(f"profile: {best_pdb.replace('.pdb', '.profile')}")
print("done")
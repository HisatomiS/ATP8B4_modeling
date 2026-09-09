import re
import numpy as np
import matplotlib.pyplot as plt

ALIGNMENT_FILE = "alignment.ali"
PROFILE_FILE = "ATP8B4_CDC50A.B99990001.profile"   # bestmodel
TEMPLATE_CODE = "8oxc"
TARGET_CODE = "ATP8B4_CDC50A"
OUT_PNG = "dope_profile_annotated.png"


def read_pir_sequences(path):
    seqs = {}
    code = None
    buf = []
    with open(path) as f:
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
    for k in seqs:
        seqs[k] = seqs[k].rstrip("*")
    return seqs


seqs = read_pir_sequences(ALIGNMENT_FILE)
templ = seqs[TEMPLATE_CODE]
targ = seqs[TARGET_CODE]

assert len(templ) == len(targ), (
    f"Length mismatch after parsing: template={len(templ)}, target={len(targ)}. "
)


gap_residues = []
target_resnum = 0
for t_char, q_char in zip(templ, targ):
    if q_char != "-":
        target_resnum += 1
        if t_char == "-":
            gap_residues.append(target_resnum)

ranges = []
if gap_residues:
    start = prev = gap_residues[0]
    for r in gap_residues[1:]:
        if r == prev + 1:
            prev = r
        else:
            ranges.append((start, prev))
            start = prev = r
    ranges.append((start, prev))

print("no residues on template:")
for s, e in ranges:
    print(f"  {s} - {e}  ({e - s + 1} residues)")

# --- plot based on DOPE data ---
data = np.loadtxt(PROFILE_FILE, comments="#", usecols=(0, -1))

plt.figure(figsize=(16, 5))
plt.plot(data[:, 0], data[:, 1], color="tab:blue", linewidth=1)

for s, e in ranges:
    plt.axvspan(s, e, color="red", alpha=0.15)

plt.xlabel("Residue index (target, ungapped)")
plt.ylabel("Normalized DOPE score")
plt.title("DOPE profile with template-uncovered (ab initio) regions shaded")
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=150)
print(f"saved {OUT_PNG}")
ALIGNMENT_FILE = "alignment.ali"
TEMPLATE_CODE = "8oxc"
TARGET_CODE = "ATP8B4_CDC50A"
TEMPLATE_START_RESNUM = 63

QUERY_ATP8B4_RESIDUES = [395, 815]


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
assert len(templ) == len(targ)

templ_resnum = TEMPLATE_START_RESNUM - 1
targ_resnum = 0
targ_to_templ = {}

for t_char, q_char in zip(templ, targ):
    if t_char != "-":
        templ_resnum += 1
    if q_char != "-":
        targ_resnum += 1
        targ_to_templ[targ_resnum] = (q_char, t_char, templ_resnum if t_char != "-" else None)

print(f"{'ATP8B4 resnum':>14} {'AA':>4}  ->  {'ATP8B1 resnum':>14} {'AA':>4}")
for resnum in QUERY_ATP8B4_RESIDUES:
    if resnum not in targ_to_templ:
        print(f"{resnum:>14}  (no mapping)")
        continue
    q_char, t_char, t_resnum = targ_to_templ[resnum]
    if t_resnum is None:
        print(f"{resnum:>14} {q_char:>4}  ->  (gap in template, no structure here)")
    else:
        print(f"{resnum:>14} {q_char:>4}  ->  {t_resnum:>14} {t_char:>4}")

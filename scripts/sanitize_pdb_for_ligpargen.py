#!/usr/bin/env python3

# Author: Star



from pathlib import Path

inp = Path("molecules/PTFMA.pdb")
out = Path("molecules/PTFMA_ligpargen.pdb")

keep = {"ATOM", "HETATM"}
serial = 1

def guess_element(atom_name: str) -> str:
    s = atom_name.strip()
    if not s:
        return "C"
    # 常见：C1, H10, O05, F0B...
    ch = s[0]
    if ch.isdigit() and len(s) > 1:
        ch = s[1]
    # 支持两字母元素（这里基本用不到）
    if len(s) >= 2 and s[:2].isalpha() and s[1].islower():
        return s[:2].title()
    return ch.upper()

lines = inp.read_text().splitlines()
out_lines = []

for ln in lines:
    if len(ln) < 6:
        continue
    rec = ln[0:6].strip()
    if rec not in keep:
        continue

    # PDB 固定列
    atom_name = ln[12:16].strip() or f"C{serial}"
    resname   = "LIG"     # 不用 UNK，改成 LIG
    chainID   = "A"
    resseq    = 1
    x = float(ln[30:38])
    y = float(ln[38:46])
    z = float(ln[46:54])
    element = (ln[76:78].strip() or guess_element(atom_name)).rjust(2)

    # 统一用 HETATM（小分子/聚合物都行）
    out_lines.append(
        f"HETATM{serial:5d} {atom_name:<4s} {resname:>3s} {chainID}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element}"
    )
    serial += 1

out_lines.append("END")
out.write_text("\n".join(out_lines) + "\n")
print("Wrote:", out)

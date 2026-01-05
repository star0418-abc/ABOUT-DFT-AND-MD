#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star




"""
make_packmol_from_recipe.py  (robust version)

Goal:
- Parse recipe YAML (supports both `recipe:` root and flat root)
- Support either:
  A) grouped format: salt_solution / polymer_matrix / crosslinker / photoinitiator
  B) flat format: components (+ optional salt dict)
- Compute target box from target_density_g_cm3 (cubic) using total mass from counts & MW
- Apply box_scale to get a looser Packmol box
- Generate Packmol input (gel.inp) + run Packmol robustly (NO stdin PIPE -> avoids Fortran 'Illegal seek')
- Emit summary.json for bookkeeping

Notes:
- If some component lacks mw_g_mol, box size cannot be computed reliably.
  Default behavior: ERROR. You can override with --allow-missing-mw (not recommended).
"""

from __future__ import annotations
import argparse
import json
import math
import os
import re
import shutil
import shlex
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml


# -----------------------------
# Small helpers
# -----------------------------
NA = 6.02214076e23  # mol^-1

CHEM_CN = {
    # 常见缩写 -> 中文全称（按你偏好：出现缩写时尽量也给中文）
    "PEGDA": "聚乙二醇二丙烯酸酯",
    "TMPTA": "三羟甲基丙烷三丙烯酸酯",
    "PTFMA": "聚(2,2,2-三氟乙基甲基丙烯酸酯)",
    "TFSI": "双(三氟甲磺酰)亚胺",
    "LiTFSI": "双(三氟甲磺酰)亚胺锂",
    "FEP": "乙基 3,3,3-三氟丙酸酯",
    "FEC": "氟代乙烯碳酸酯",
    "EC": "碳酸乙烯酯",
}

def eprint(*args):
    print(*args, file=sys.stderr)

def slugify(s: str) -> str:
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s or "item"

def read_yaml(path: Path) -> Dict[str, Any]:
    try:
        txt = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        # 兜底
        txt = path.read_text(errors="replace")
    data = yaml.safe_load(txt)
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping/dict, got: {type(data)}")
    return data

def find_packmol_exe() -> Optional[str]:
    exe = shutil.which("packmol")
    if exe:
        return exe
    # common fallbacks
    for p in ["/usr/bin/packmol", "/usr/local/bin/packmol", str(Path.home() / ".local/bin/packmol")]:
        if Path(p).exists():
            return p
    return None

def resolve_file(path_str: str, cfg_path: Path, project_root: Path) -> Path:
    """
    Resolve a file path appearing in YAML.

    Search order:
    1) absolute path
    2) relative to config file dir
    3) relative to project root
    4) project_root/molecules/<basename>
    """
    p = Path(path_str)
    if p.is_absolute() and p.exists():
        return p

    cand = (cfg_path.parent / p)
    if cand.exists():
        return cand

    cand = (project_root / p)
    if cand.exists():
        return cand

    cand = (project_root / "molecules" / p.name)
    if cand.exists():
        return cand

    raise FileNotFoundError(f"Cannot find file '{path_str}' (tried cfg_dir, project_root, molecules/)")

def safe_float(x: Any, default: float) -> float:
    try:
        return float(x)
    except Exception:
        return default


# -----------------------------
# Data models
# -----------------------------
@dataclass
class Structure:
    name: str
    pdb_src: Path          # original source path
    pdb_local: Path        # copied-to packmol dir path
    count: int
    mw_g_mol: float
    charge: float = 0.0
    kind: str = "component"  # component / salt_cation / salt_anion

@dataclass
class BoxInfo:
    shape: str
    Lx_A: float
    Ly_A: float
    Lz_A: float

    @property
    def volume_A3(self) -> float:
        return self.Lx_A * self.Ly_A * self.Lz_A

    @property
    def volume_nm3(self) -> float:
        return self.volume_A3 / 1000.0

    @property
    def volume_cm3(self) -> float:
        return self.volume_A3 * 1e-24  # 1 A^3 = 1e-24 cm^3

    @property
    def volume_L(self) -> float:
        return self.volume_cm3 * 1e-3  # 1 L = 1000 cm^3


# -----------------------------
# Parsing recipe
# -----------------------------
def get_root_cfg(data: Dict[str, Any]) -> Dict[str, Any]:
    # allow either top-level recipe or flat
    root = data.get("recipe", data)
    if not isinstance(root, dict):
        raise ValueError("recipe content must be a dict/mapping")
    return root

def parse_salt_from_salt_dict(salt: Dict[str, Any], cfg_path: Path, project_root: Path) -> Tuple[List[Structure], Dict[str, Any]]:
    """
    Parse a `salt:` dict of the style:
      salt:
        name: LiTFSI
        cation: {file: ..., mw_g_mol: 6.941, charge: +1, stoichiometry: 1, count: 474}
        anion:  {file: ..., mw_g_mol: 280.15, charge: -1, stoichiometry: 1, count: 474}
        n_formula_units: 474   (optional)
        count: 474            (optional)
    """
    salt_name = str(salt.get("name", "salt"))
    n_formula = salt.get("n_formula_units", salt.get("count", None))
    cation = salt.get("cation", {}) or {}
    anion = salt.get("anion", {}) or {}

    # infer n_formula from cation/anion counts if needed
    if n_formula is None:
        # try: min(cation.count/stoich, anion.count/stoich)
        try:
            c_cnt = int(cation.get("count", 0))
            a_cnt = int(anion.get("count", 0))
            c_sto = int(cation.get("stoichiometry", 1))
            a_sto = int(anion.get("stoichiometry", 1))
            n_formula = min(c_cnt // max(c_sto, 1), a_cnt // max(a_sto, 1))
        except Exception:
            n_formula = 0

    n_formula = int(n_formula or 0)
    if n_formula <= 0:
        raise ValueError(f"Salt '{salt_name}' has no valid n_formula_units/count (>0).")

    def make_ion(ion: Dict[str, Any], kind: str) -> Structure:
        ion_name = str(ion.get("name", f"{salt_name}_{kind}"))
        ion_file = resolve_file(str(ion.get("file", "")), cfg_path, project_root)
        mw = safe_float(ion.get("mw_g_mol", 0.0), 0.0)
        charge = safe_float(ion.get("charge", 0.0), 0.0)
        sto = int(ion.get("stoichiometry", 1) or 1)
        cnt = n_formula * sto
        return Structure(
            name=ion_name,
            pdb_src=ion_file,
            pdb_local=Path(ion_file.name),
            count=cnt,
            mw_g_mol=mw,
            charge=charge,
            kind=kind,
        )

    structures = [make_ion(cation, "salt_cation"), make_ion(anion, "salt_anion")]
    meta = {
        "name": salt_name,
        "n_formula_units": n_formula,
        "cation": {"name": structures[0].name, "count": structures[0].count, "mw_g_mol": structures[0].mw_g_mol, "charge": structures[0].charge},
        "anion":  {"name": structures[1].name, "count": structures[1].count, "mw_g_mol": structures[1].mw_g_mol, "charge": structures[1].charge},
    }
    return structures, meta

def parse_grouped(root: Dict[str, Any], cfg_path: Path, project_root: Path) -> Tuple[List[Structure], Dict[str, Any]]:
    """
    Parse grouped format: salt_solution / polymer_matrix / crosslinker / photoinitiator
    Return: structures list + salt meta (optional)
    """
    structures: List[Structure] = []
    salt_meta: Dict[str, Any] = {}

    # --- salt_solution first (按你偏好顺序：盐溶液->聚合物基质->交联剂->引发剂)
    if "salt_solution" in root and isinstance(root["salt_solution"], list):
        for item in root["salt_solution"]:
            if not isinstance(item, dict):
                continue
            itype = str(item.get("type", "solvent"))
            if itype == "salt":
                # formula units count
                n_formula = int(item.get("count", 0) or 0)
                if n_formula <= 0:
                    raise ValueError(f"salt_solution salt '{item.get('name','salt')}' count must be > 0")

                salt_name = str(item.get("name", "salt"))

                cation = item.get("cation", {}) or {}
                anion  = item.get("anion", {}) or {}

                # build salt dict and reuse the other parser
                salt_dict = {
                    "name": salt_name,
                    "count": n_formula,
                    "cation": cation,
                    "anion": anion,
                }
                salt_structs, meta = parse_salt_from_salt_dict(salt_dict, cfg_path, project_root)
                structures.extend(salt_structs)
                salt_meta = meta
            else:
                name = str(item.get("name", item.get("id", "solvent")))
                pdb = resolve_file(str(item.get("file", "")), cfg_path, project_root)
                cnt = int(item.get("count", 0) or 0)
                mw = safe_float(item.get("mw_g_mol", 0.0), 0.0)
                charge = safe_float(item.get("charge", 0.0), 0.0)
                if cnt <= 0:
                    continue
                structures.append(
                    Structure(
                        name=name, pdb_src=pdb, pdb_local=Path(pdb.name),
                        count=cnt, mw_g_mol=mw, charge=charge, kind="component"
                    )
                )

    # --- polymer_matrix
    if "polymer_matrix" in root and isinstance(root["polymer_matrix"], list):
        for item in root["polymer_matrix"]:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name", item.get("id", "polymer")))
            pdb = resolve_file(str(item.get("file", "")), cfg_path, project_root)
            cnt = int(item.get("count", 0) or 0)
            mw = safe_float(item.get("mw_g_mol", 0.0), 0.0)
            charge = safe_float(item.get("charge", 0.0), 0.0)
            if cnt <= 0:
                continue
            structures.append(
                Structure(name=name, pdb_src=pdb, pdb_local=Path(pdb.name),
                          count=cnt, mw_g_mol=mw, charge=charge, kind="component")
            )

    # --- crosslinker (optional)
    if "crosslinker" in root and isinstance(root["crosslinker"], list):
        for item in root["crosslinker"]:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name", item.get("id", "crosslinker")))
            pdb = resolve_file(str(item.get("file", "")), cfg_path, project_root)
            cnt = int(item.get("count", 0) or 0)
            mw = safe_float(item.get("mw_g_mol", 0.0), 0.0)
            charge = safe_float(item.get("charge", 0.0), 0.0)
            if cnt <= 0:
                continue
            structures.append(
                Structure(name=name, pdb_src=pdb, pdb_local=Path(pdb.name),
                          count=cnt, mw_g_mol=mw, charge=charge, kind="component")
            )

    # --- photoinitiator (optional)
    if "photoinitiator" in root and isinstance(root["photoinitiator"], list):
        for item in root["photoinitiator"]:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name", item.get("id", "initiator")))
            pdb = resolve_file(str(item.get("file", "")), cfg_path, project_root)
            cnt = int(item.get("count", 0) or 0)
            mw = safe_float(item.get("mw_g_mol", 0.0), 0.0)
            charge = safe_float(item.get("charge", 0.0), 0.0)
            if cnt <= 0:
                continue
            structures.append(
                Structure(name=name, pdb_src=pdb, pdb_local=Path(pdb.name),
                          count=cnt, mw_g_mol=mw, charge=charge, kind="component")
            )

    return structures, salt_meta

def parse_flat(root: Dict[str, Any], cfg_path: Path, project_root: Path) -> Tuple[List[Structure], Dict[str, Any]]:
    """
    Parse flat format:
      - components: [ {name, file, count, mw_g_mol, charge, type}, ...]
      - optional salt: {...}  (will be decomposed)
    """
    structures: List[Structure] = []
    salt_meta: Dict[str, Any] = {}

    # optional salt dict
    if "salt" in root and isinstance(root["salt"], dict):
        salt_structs, meta = parse_salt_from_salt_dict(root["salt"], cfg_path, project_root)
        structures.extend(salt_structs)
        salt_meta = meta

    # components list
    comps = root.get("components", None)
    if not isinstance(comps, list):
        raise ValueError("Flat format requires 'components' as a list.")

    for item in comps:
        if not isinstance(item, dict):
            continue
        itype = str(item.get("type", "component"))
        if itype == "salt":
            # allow salt inside components too
            salt_structs, meta = parse_salt_from_salt_dict(item, cfg_path, project_root)
            structures.extend(salt_structs)
            salt_meta = meta
            continue

        name = str(item.get("name", item.get("id", "component")))
        pdb = resolve_file(str(item.get("file", "")), cfg_path, project_root)
        cnt = int(item.get("count", 0) or 0)
        mw = safe_float(item.get("mw_g_mol", 0.0), 0.0)
        charge = safe_float(item.get("charge", 0.0), 0.0)
        if cnt <= 0:
            continue
        structures.append(
            Structure(name=name, pdb_src=pdb, pdb_local=Path(pdb.name),
                      count=cnt, mw_g_mol=mw, charge=charge, kind="component")
        )

    return structures, salt_meta

def parse_recipe(cfg_path: Path, allow_missing_mw: bool) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], List[Structure], Dict[str, Any]]:
    """
    Return: (root, system_cfg, packmol_cfg, structures, salt_meta)
    """
    data = read_yaml(cfg_path)
    root = get_root_cfg(data)

    system_cfg = root.get("system", {}) or {}
    packmol_cfg = root.get("packmol", {}) or {}
    constraints = root.get("constraints", {}) or root.get("constraint", {}) or {}

    # defaults
    packmol_cfg.setdefault("tolerance_A", 2.0)
    packmol_cfg.setdefault("box_scale", 1.5)
    packmol_cfg.setdefault("seed", 12345)
    packmol_cfg.setdefault("filetype", "pdb")
    packmol_cfg.setdefault("output_pdb", "gel.pdb")

    # parse structures
    has_grouped = any(k in root for k in ("salt_solution", "polymer_matrix", "crosslinker", "photoinitiator"))
    structures: List[Structure] = []
    salt_meta: Dict[str, Any] = {}

    if has_grouped:
        structures, salt_meta = parse_grouped(root, cfg_path, project_root=cfg_path.parent.parent)
    else:
        structures, salt_meta = parse_flat(root, cfg_path, project_root=cfg_path.parent.parent)

    if not structures:
        raise ValueError("No structures parsed (all counts are zero?).")

    # MW check
    if not allow_missing_mw:
        missing = [s for s in structures if (s.mw_g_mol is None or s.mw_g_mol <= 0)]
        if missing:
            msg = "Missing/invalid mw_g_mol for:\n" + "\n".join([f"  - {m.name} ({m.pdb_src})" for m in missing])
            msg += "\nFix recipe YAML by adding mw_g_mol, or re-run with --allow-missing-mw (not recommended)."
            raise ValueError(msg)

    return root, system_cfg, packmol_cfg, constraints, structures, salt_meta


# -----------------------------
# Box & Packmol I/O
# -----------------------------
def compute_target_box_from_density(structures: List[Structure], target_density_g_cm3: float) -> BoxInfo:
    total_mass_g_mol = sum(s.count * float(s.mw_g_mol or 0.0) for s in structures)
    total_mass_g = total_mass_g_mol / NA
    volume_cm3 = total_mass_g / target_density_g_cm3
    L_cm = volume_cm3 ** (1.0 / 3.0)
    L_A = L_cm * 1e8
    return BoxInfo(shape="cubic", Lx_A=L_A, Ly_A=L_A, Lz_A=L_A)

def compute_packmol_box(target: BoxInfo, box_scale: float) -> BoxInfo:
    return BoxInfo(shape=target.shape,
                   Lx_A=target.Lx_A * box_scale,
                   Ly_A=target.Ly_A * box_scale,
                   Lz_A=target.Lz_A * box_scale)

def write_packmol_inp(inp_path: Path, output_pdb: str, tolerance_A: float, seed: int, box: BoxInfo, structures: List[Structure]) -> None:
    lines: List[str] = []
    lines.append(f"tolerance {tolerance_A}")
    lines.append("filetype pdb")
    lines.append(f"output {output_pdb}")
    lines.append(f"seed {int(seed)}")
    lines.append("")  # blank

    Lx, Ly, Lz = box.Lx_A, box.Ly_A, box.Lz_A

    for s in structures:
        # use local file name in packmol dir
        lines.append(f"structure {s.pdb_local.name}")
        lines.append(f"  number {int(s.count)}")
        lines.append(f"  inside box 0. 0. 0. {Lx:.4f} {Ly:.4f} {Lz:.4f}")
        lines.append("end structure")
        lines.append("")

    inp_path.write_text("\n".join(lines), encoding="utf-8")

def copy_structures_to_dir(structures: List[Structure], out_dir: Path) -> None:
    """
    Copy all pdb_src to out_dir and set pdb_local accordingly (unique names if needed).
    """
    used = set()
    for s in structures:
        base = slugify(s.pdb_src.stem)
        ext = s.pdb_src.suffix or ".pdb"
        # name collisions: add component name prefix if needed
        preferred = slugify(s.name)
        fname = f"{preferred}{ext}"
        if fname in used:
            fname = f"{preferred}_{base}{ext}"
        if fname in used:
            # fallback with numeric suffix
            k = 2
            while f"{preferred}_{k}{ext}" in used:
                k += 1
            fname = f"{preferred}_{k}{ext}"

        dst = out_dir / fname
        shutil.copy2(s.pdb_src, dst)
        used.add(fname)
        s.pdb_local = dst.name and Path(dst.name)  # store local relative name

def run_packmol(packmol_exe: str, inp_path: Path, log_path: Path, cwd: Path, output_pdb: str = "gel.pdb") -> None:
    """
    使用 shell 重定向方式运行 Packmol，避免 Fortran 的 'Illegal seek' 错误。
    
    关键：不使用 stdin=PIPE 或 input=...，而是用 shell 的 < 重定向，
    让 Fortran 的 stdin 来自真实文件（可 seek）。
    """
    # 使用 shlex.quote 转义路径，避免空格问题
    cmd_str = (
        f"{shlex.quote(packmol_exe)} "
        f"< {shlex.quote(str(inp_path))} "
        f"> {shlex.quote(str(log_path))} 2>&1"
    )
    print(f"[CMD] {cmd_str}")
    
    ret = subprocess.call(
        cmd_str,
        shell=True,
        cwd=str(cwd),
        executable="/bin/bash",
    )
    
    if ret != 0:
        raise RuntimeError(f"Packmol failed (code={ret}). Check log: {log_path}")
    
    # 检查输出文件
    out_pdb = cwd / output_pdb
    if not out_pdb.exists():
        raise RuntimeError(f"Packmol finished but output PDB not found: {out_pdb}. Check log: {log_path}")
    if out_pdb.stat().st_size < 2000:
        raise RuntimeError(f"Packmol finished but output PDB abnormally small ({out_pdb.stat().st_size} bytes). Check log: {log_path}")

def compute_mass_summary(structures: List[Structure]) -> Dict[str, Any]:
    total_mass_g_mol = sum(s.count * float(s.mw_g_mol or 0.0) for s in structures)
    total_mass_g = total_mass_g_mol / NA
    total_mol = sum(s.count for s in structures) / NA  # not chemically meaningful for mixture, but ok as a rough "molecule count mol"
    total_molecules = int(sum(s.count for s in structures))
    return {
        "total_mass_g": total_mass_g,
        "total_mass_g_mol_equiv": total_mass_g_mol,
        "total_molecules": total_molecules,
        "note": "total_mass_g computed from sum(count*mw)/NA; total_molecules is total number of structures placed in Packmol",
    }

def compute_charge(structures: List[Structure]) -> Dict[str, Any]:
    net = sum(float(s.charge or 0.0) * int(s.count) for s in structures)
    return {"net_charge_e": net}

def write_summary_json(
    path: Path,
    recipe_name: str,
    system_cfg: Dict[str, Any],
    packmol_cfg: Dict[str, Any],
    constraints: Dict[str, Any],
    target_box: BoxInfo,
    packmol_box: BoxInfo,
    structures: List[Structure],
    salt_meta: Dict[str, Any],
    output_pdb: str,
) -> None:
    mass = compute_mass_summary(structures)
    charge = compute_charge(structures)

    # packmol initial density based on packmol box
    total_mass_g = float(mass["total_mass_g"])
    rho_packmol = total_mass_g / packmol_box.volume_cm3 if packmol_box.volume_cm3 > 0 else None

    components_out = []
    for s in structures:
        components_out.append({
            "name": s.name,
            "file": str(s.pdb_local),
            "count": int(s.count),
            "mw_g_mol": float(s.mw_g_mol or 0.0),
            "charge": float(s.charge or 0.0),
            "kind": s.kind,
        })

    out = {
        "recipe_name": recipe_name,
        "generator_version": "2.2",
        "system": system_cfg,
        "input_parameters": {
            "target_density_g_cm3": safe_float(system_cfg.get("target_density_g_cm3", 1.0), 1.0),
            "box_shape": system_cfg.get("box_shape", "cubic"),
            "tolerance_A": safe_float(packmol_cfg.get("tolerance_A", 2.0), 2.0),
            "packmol_seed": int(packmol_cfg.get("seed", 12345)),
            "packmol_box_scale": safe_float(packmol_cfg.get("box_scale", 1.5), 1.5),
        },
        "target_box": {
            "shape": target_box.shape,
            "Lx_A": target_box.Lx_A, "Ly_A": target_box.Ly_A, "Lz_A": target_box.Lz_A,
            "volume_A3": target_box.volume_A3,
            "volume_nm3": target_box.volume_nm3,
            "volume_cm3": target_box.volume_cm3,
            "volume_L": target_box.volume_L,
        },
        "packmol_box": {
            "Lx_A": packmol_box.Lx_A, "Ly_A": packmol_box.Ly_A, "Lz_A": packmol_box.Lz_A,
            "volume_A3": packmol_box.volume_A3,
            "volume_nm3": packmol_box.volume_nm3,
            "volume_cm3": packmol_box.volume_cm3,
            "initial_density_g_cm3": rho_packmol,
            "note": "Packmol box is intentionally looser; target density is achieved later by GROMACS NPT.",
        },
        "mass_summary": mass,
        "charge_summary": charge,
        "constraints_applied": constraints,
        "salt": salt_meta or None,
        "structures": components_out,
        "output_files": {
            "packmol_input": "gel.inp",
            "packmol_output": output_pdb,
        },
        "abbr_cn": CHEM_CN,  # 方便你查缩写中文
    }

    path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Generate & run Packmol from recipe YAML (robust).")
    ap.add_argument("-c", "--config", required=True, help="Recipe YAML path (e.g., configs/recipe_wsgpe_repro.yaml)")
    ap.add_argument("-o", "--out", required=True, help="Output directory (e.g., outputs/WSGPE_Reproduction_Final)")
    ap.add_argument("--dry-run", action="store_true", help="Only generate gel.inp + summary.json, do not run packmol")
    ap.add_argument("--keep", action="store_true", help="Do not delete existing output directory")
    ap.add_argument("--allow-missing-mw", action="store_true", help="Allow missing mw_g_mol (box size may be wrong)")
    ap.add_argument("--tolerance", type=float, default=None, help="Override packmol tolerance (A)")
    ap.add_argument("--box-scale", type=float, default=None, help="Override packmol box_scale")
    ap.add_argument("--seed", type=int, default=None, help="Override packmol seed")
    ap.add_argument("--output-pdb", type=str, default=None, help="Override output pdb name (default: gel.pdb)")
    args = ap.parse_args()

    cfg_path = Path(args.config).resolve()
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")

    out_root = Path(args.out).resolve()
    packmol_dir = out_root / "packmol"

    if out_root.exists() and not args.keep:
        shutil.rmtree(out_root)
    packmol_dir.mkdir(parents=True, exist_ok=True)

    # parse
    root, system_cfg, packmol_cfg, constraints, structures, salt_meta = parse_recipe(cfg_path, allow_missing_mw=args.allow_missing_mw)

    # overrides
    if args.tolerance is not None:
        packmol_cfg["tolerance_A"] = float(args.tolerance)
    if args.box_scale is not None:
        packmol_cfg["box_scale"] = float(args.box_scale)
    if args.seed is not None:
        packmol_cfg["seed"] = int(args.seed)
    if args.output_pdb is not None:
        packmol_cfg["output_pdb"] = str(args.output_pdb)

    recipe_name = str(system_cfg.get("name", out_root.name))
    target_density = safe_float(system_cfg.get("target_density_g_cm3", 1.0), 1.0)
    box_shape = str(system_cfg.get("box_shape", "cubic"))

    if box_shape.lower() != "cubic":
        raise NotImplementedError("Only cubic box is supported in this script for now.")

    # copy pdbs into packmol dir for robust relative paths
    copy_structures_to_dir(structures, packmol_dir)

    # compute boxes
    target_box = compute_target_box_from_density(structures, target_density_g_cm3=target_density)
    box_scale = safe_float(packmol_cfg.get("box_scale", 1.5), 1.5)
    pack_box = compute_packmol_box(target_box, box_scale=box_scale)

    # generate input
    inp_path = packmol_dir / "gel.inp"
    output_pdb = str(packmol_cfg.get("output_pdb", "gel.pdb"))
    tol = safe_float(packmol_cfg.get("tolerance_A", 2.0), 2.0)
    seed = int(packmol_cfg.get("seed", 12345))

    write_packmol_inp(inp_path, output_pdb=output_pdb, tolerance_A=tol, seed=seed, box=pack_box, structures=structures)

    # summary
    summary_path = packmol_dir / "summary.json"
    write_summary_json(
        summary_path,
        recipe_name=recipe_name,
        system_cfg=system_cfg,
        packmol_cfg=packmol_cfg,
        constraints=constraints,
        target_box=target_box,
        packmol_box=pack_box,
        structures=structures,
        salt_meta=salt_meta,
        output_pdb=output_pdb,
    )

    print(f"[INFO] Config: {cfg_path}")
    print(f"[INFO] Output: {out_root}")
    print(f"[INFO] Packmol dir: {packmol_dir}")
    print(f"[INFO] gel.inp: {inp_path}")
    print(f"[INFO] summary.json: {summary_path}")
    print(f"[INFO] target box (A): {target_box.Lx_A:.4f}")
    print(f"[INFO] packmol box (A): {pack_box.Lx_A:.4f} (scale={box_scale})")

    if args.dry_run:
        print("[INFO] dry-run enabled: skipping Packmol execution.")
        return

    packmol_exe = find_packmol_exe()
    if not packmol_exe:
        raise RuntimeError("Cannot find packmol executable. Install it or add to PATH.")

    log_path = packmol_dir / "packmol.log"
    print(f"[RUN] {packmol_exe} < gel.inp  (cwd={packmol_dir})")
    run_packmol(packmol_exe, inp_path=inp_path, log_path=log_path, cwd=packmol_dir, output_pdb=output_pdb)
    
    out_pdb = packmol_dir / output_pdb

    print("[SUCCESS] Packmol done.")
    print(f"  -> {out_pdb} ({out_pdb.stat().st_size/1024/1024:.2f} MB)")
    print("Next: run your GROMACS stage (run_gmx.sh / run_full.sh) using this output directory.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        eprint(f"[FATAL] {e}")
        sys.exit(1)

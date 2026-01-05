#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
generate_topology_gaff.py - 使用 GAFF2 力场生成 GROMACS 拓扑
================================================================
为电解质体系生成完整的 GROMACS 拓扑，使用通用 AMBER 力场 (GAFF2)。

策略:
1. 使用 acpype 为每个分子生成 GAFF2 参数
2. 合并所有 atomtypes 到一个文件
3. 创建主 topol.top

用法:
    python3 scripts/generate_topology_gaff.py \
        -i outputs/WSGPE_Reproduction_Final/packmol \
        -c configs/recipe_wsgpe_repro.yaml

"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml


def parse_recipe(recipe_path: Path) -> dict:
    """解析配方文件"""
    with open(recipe_path) as f:
        cfg = yaml.safe_load(f)
    
    recipe = cfg.get('recipe', cfg)
    components = recipe.get('components', [])
    
    molecules = {}
    for comp in components:
        mol_id = comp['id']
        molecules[mol_id] = {
            'name': comp['name'],
            'file': comp['file'],
            'count': comp['count'],
            'charge': int(comp.get('charge', 0)),
        }
    
    return molecules


def run_acpype_for_mol(mol_pdb: Path, mol_name: str, charge: int, work_dir: Path) -> Optional[Path]:
    """为单个分子运行 acpype"""
    print(f"  [acpype] {mol_name} (charge={charge})...")
    
    acpype_dir = work_dir / f"acpype_{mol_name}"
    acpype_dir.mkdir(parents=True, exist_ok=True)
    
    # 复制 PDB
    pdb_copy = acpype_dir / f"{mol_name}.pdb"
    shutil.copy(mol_pdb.resolve(), pdb_copy)
    
    # 运行 acpype
    cmd = [
        "acpype",
        "-i", str(pdb_copy.resolve()),
        "-n", str(charge),
        "-a", "gaff2",
        "-o", "gmx",
        "-b", mol_name,
    ]
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(acpype_dir),
            capture_output=True,
            text=True,
            timeout=600,
        )
    except Exception as e:
        print(f"    [ERROR] {e}")
        return None
    
    # 查找 itp 文件
    acpype_out = acpype_dir / f"{mol_name}.acpype"
    itp_file = acpype_out / f"{mol_name}_GMX.itp"
    
    if itp_file.exists():
        return itp_file
    
    # 尝试其他位置
    for itp in acpype_dir.rglob("*_GMX.itp"):
        return itp
    for itp in acpype_dir.rglob("*.itp"):
        return itp
    
    print(f"    [WARN] 未找到 itp: {result.stdout[-200:]}")
    return None


def extract_atomtypes(itp_content: str) -> str:
    """从 itp 内容中提取 atomtypes 部分"""
    lines = itp_content.split('\n')
    atomtypes = []
    in_atomtypes = False
    
    for line in lines:
        if '[ atomtypes ]' in line:
            in_atomtypes = True
            continue
        elif in_atomtypes and line.strip().startswith('['):
            break
        elif in_atomtypes and line.strip() and not line.strip().startswith(';'):
            atomtypes.append(line)
    
    return '\n'.join(atomtypes)


def extract_moleculetype(itp_content: str, new_name: str) -> str:
    """提取 moleculetype 及后续部分，跳过 atomtypes"""
    lines = itp_content.split('\n')
    result = []
    skip_atomtypes = False
    started = False
    
    for line in lines:
        if '[ atomtypes ]' in line:
            skip_atomtypes = True
            continue
        elif skip_atomtypes and line.strip().startswith('['):
            skip_atomtypes = False
        
        if skip_atomtypes:
            continue
        
        if '[ moleculetype ]' in line:
            started = True
        
        if started:
            result.append(line)
    
    return '\n'.join(result)


def create_li_itp(output_dir: Path) -> Tuple[str, str]:
    """创建 Li+ 离子的 itp（使用标准 GAFF2 参数）"""
    atomtypes = """Li      Li          0.00000  6.94000   A   1.2500e-01  6.2760e-02"""
    
    moleculetype = """; Li+ ion topology

[ moleculetype ]
; name  nrexcl
Li  3

[ atoms ]
;   nr  type  resnr  residue  atom  cgnr  charge   mass
     1  Li       1    LI+     LI     1   1.0000   6.94
"""
    return atomtypes, moleculetype


def create_simple_itp(pdb_file: Path, mol_name: str, charge: int) -> Tuple[str, str]:
    """创建简化的 itp（基于 GAFF2 原子类型）"""
    # 读取 PDB
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                atoms.append((atom_name, element))
    
    if not atoms:
        return "", ""
    
    # GAFF2 原子类型映射
    gaff_map = {
        'C': 'c3',   # sp3 carbon
        'H': 'hc',   # H on carbon
        'O': 'o',    # carbonyl O
        'N': 'n3',   # sp3 nitrogen
        'S': 's4',   # sulfur
        'F': 'f',    # fluorine
        'P': 'p5',   # phosphorus
        'Li': 'Li',  # lithium
    }
    
    # GAFF2 原子类型参数 (sigma, epsilon)
    gaff_params = {
        'c3': (3.39967e-01, 4.57730e-01),
        'hc': (2.64953e-01, 6.56888e-02),
        'o': (2.95992e-01, 8.78640e-01),
        'os': (3.00001e-01, 7.11280e-01),
        'n3': (3.25000e-01, 7.11280e-01),
        's4': (3.56359e-01, 1.04600e+00),
        'f': (3.11815e-01, 2.55224e-01),
        'p5': (3.74177e-01, 8.36800e-01),
    }
    
    # 原子质量
    mass_map = {
        'C': 12.01, 'H': 1.008, 'O': 16.00, 'N': 14.01,
        'S': 32.06, 'F': 19.00, 'P': 30.97, 'Li': 6.94,
    }
    
    # 生成 atomtypes
    atomtypes_set = set()
    atomtypes_lines = []
    
    for atom_name, element in atoms:
        atype = gaff_map.get(element, 'c3')
        if atype not in atomtypes_set:
            atomtypes_set.add(atype)
            params = gaff_params.get(atype, (3.4e-01, 4.5e-01))
            mass = mass_map.get(element, 12.0)
            atomtypes_lines.append(
                f"{atype:4s}   {atype:4s}  0.00000  {mass:.3f}   A   {params[0]:.5e}  {params[1]:.5e}"
            )
    
    atomtypes = '\n'.join(atomtypes_lines)
    
    # 生成 moleculetype
    charge_per_atom = charge / len(atoms) if atoms else 0
    
    mol_lines = [
        f"; {mol_name} topology",
        "",
        "[ moleculetype ]",
        "; name  nrexcl",
        f"{mol_name}  3",
        "",
        "[ atoms ]",
        ";   nr  type  resnr  residue  atom  cgnr  charge   mass",
    ]
    
    for i, (atom_name, element) in enumerate(atoms, 1):
        atype = gaff_map.get(element, 'c3')
        mass = mass_map.get(element, 12.0)
        mol_lines.append(
            f"     {i:3d}  {atype:4s}    1  {mol_name:4s}  {atom_name:4s}  {i:3d}  {charge_per_atom:8.4f}  {mass:.2f}"
        )
    
    # 添加简单的键（根据距离或连接信息）
    # 这里暂时不添加键信息，让 GROMACS 使用默认处理
    
    moleculetype = '\n'.join(mol_lines)
    
    return atomtypes, moleculetype


def create_topol_top(
    molecules: dict,
    mol_data: Dict[str, Tuple[str, str]],  # {mol_name: (atomtypes, moleculetype)}
    output_dir: Path,
) -> Path:
    """创建主 topol.top 文件"""
    topol_path = output_dir / "topol.top"
    
    # 合并所有 atomtypes
    all_atomtypes = set()
    atomtypes_lines = []
    
    for mol_name, (atomtypes, _) in mol_data.items():
        for line in atomtypes.split('\n'):
            if line.strip() and not line.startswith(';'):
                atype = line.split()[0]
                if atype not in all_atomtypes:
                    all_atomtypes.add(atype)
                    atomtypes_lines.append(line)
    
    with open(topol_path, 'w') as f:
        f.write("; WSGPE Electrolyte Topology\n")
        f.write("; Generated by generate_topology_gaff.py\n")
        f.write("; Force field: GAFF2 (General AMBER Force Field 2)\n\n")
        
        # 默认参数
        f.write("[ defaults ]\n")
        f.write("; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
        f.write("1         2          yes        0.5      0.8333\n\n")
        
        # Atomtypes
        f.write("[ atomtypes ]\n")
        f.write("; name  bond_type  mass  charge  ptype  sigma  epsilon\n")
        for line in atomtypes_lines:
            f.write(f"{line}\n")
        f.write("\n")
        
        # Moleculetypes
        for mol_name, (_, moleculetype) in mol_data.items():
            f.write(f"\n; ========== {mol_name} ==========\n")
            f.write(moleculetype)
            f.write("\n")
        
        # System
        f.write("\n[ system ]\n")
        f.write("WSGPE Electrolyte\n\n")
        
        # Molecules
        f.write("[ molecules ]\n")
        f.write("; Compound        nmols\n")
        for mol_id, info in molecules.items():
            mol_name = info['name']
            count = info['count']
            f.write(f"{mol_name:<16s} {count}\n")
    
    return topol_path


def pdb_to_gro(pdb_file: Path, gro_file: Path, box_nm: float) -> bool:
    """PDB -> GRO"""
    cmd = [
        "gmx", "editconf",
        "-f", str(pdb_file),
        "-o", str(gro_file),
        "-box", str(box_nm), str(box_nm), str(box_nm),
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return gro_file.exists()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-o", "--output")
    parser.add_argument("--box", type=float, default=25.4)
    parser.add_argument("--skip-acpype", action="store_true", help="跳过 acpype，使用简化参数")
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    config_path = Path(args.config)
    output_dir = Path(args.output) if args.output else input_dir
    base_dir = config_path.parent.parent
    
    gel_pdb = input_dir / "gel.pdb"
    if not gel_pdb.exists():
        print(f"[ERROR] {gel_pdb} 不存在")
        return 1
    
    print(f"[INFO] 输入: {gel_pdb}")
    print(f"[INFO] 输出: {output_dir}")
    
    molecules = parse_recipe(config_path)
    print(f"[INFO] 分子数: {len(molecules)}")
    
    mol_data = {}
    
    for mol_id, info in molecules.items():
        mol_name = info['name']
        mol_file = base_dir / info['file']
        charge = info['charge']
        
        print(f"\n处理: {mol_name}")
        
        if mol_name.upper() == "LI":
            # Li+ 使用预定义参数
            atomtypes, moleculetype = create_li_itp(output_dir)
            mol_data[mol_name] = (atomtypes, moleculetype)
            print(f"  [OK] 使用预定义 Li+ 参数")
            continue
        
        if not args.skip_acpype:
            # 尝试 acpype
            itp = run_acpype_for_mol(mol_file, mol_name, charge, output_dir)
            
            if itp and itp.exists():
                content = itp.read_text()
                atomtypes = extract_atomtypes(content)
                moleculetype = extract_moleculetype(content, mol_name)
                mol_data[mol_name] = (atomtypes, moleculetype)
                print(f"  [OK] acpype: {itp}")
                continue
        
        # Fallback: 简化参数
        print(f"  [FALLBACK] 使用简化 GAFF2 参数")
        atomtypes, moleculetype = create_simple_itp(mol_file, mol_name, charge)
        mol_data[mol_name] = (atomtypes, moleculetype)
    
    # 创建 topol.top
    print(f"\n[INFO] 创建 topol.top...")
    topol = create_topol_top(molecules, mol_data, output_dir)
    print(f"[OK] {topol}")
    
    # PDB -> GRO
    print(f"[INFO] 转换 PDB -> GRO (box={args.box} nm)...")
    gro = output_dir / "system.gro"
    if pdb_to_gro(gel_pdb, gro, args.box):
        print(f"[OK] {gro}")
    
    print(f"\n[SUCCESS] 拓扑生成完成!")
    print(f"\n下一步:")
    print(f"  ./scripts/run_gmx.sh -i {output_dir} -o {output_dir.parent}/gmx")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


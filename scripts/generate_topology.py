#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
generate_topology.py - 自动生成 GROMACS 拓扑文件
===================================================
使用 acpype 为每个分子生成 GAFF 力场参数，并创建主 topol.top 文件。

用法:
    python3 scripts/generate_topology.py -i outputs/WSGPE_Reproduction_Final/packmol \
        -c configs/recipe_wsgpe_repro.yaml

"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def parse_recipe(recipe_path: Path) -> dict:
    """解析配方文件，返回分子列表"""
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


def run_acpype(pdb_file: Path, output_dir: Path, mol_name: str, charge: int = 0) -> Path:
    """运行 acpype 生成 itp 文件"""
    print(f"  [acpype] 处理 {mol_name} (charge={charge})...")
    
    # 创建工作目录（使用绝对路径）
    work_dir = Path(output_dir).resolve() / f"acpype_{mol_name}"
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # 复制 PDB 到工作目录（使用绝对路径）
    pdb_copy = work_dir / f"{mol_name}.pdb"
    shutil.copy(Path(pdb_file).resolve(), pdb_copy)
    
    # 确保 PDB 文件存在
    if not pdb_copy.exists():
        print(f"  [ERROR] PDB 复制失败: {pdb_copy}")
        return None
    
    # 运行 acpype（使用绝对路径）
    cmd = [
        "acpype",
        "-i", str(pdb_copy.resolve()),
        "-n", str(charge),
        "-a", "gaff2",  # 使用 GAFF2 力场
        "-o", "gmx",    # 输出 GROMACS 格式
        "-b", mol_name,
    ]
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=300,  # 5分钟超时
        )
        
        if result.returncode != 0:
            print(f"  [WARN] acpype 返回非零: {result.returncode}")
            print(f"         stdout: {result.stdout[-500:] if result.stdout else 'N/A'}")
            print(f"         stderr: {result.stderr[-500:] if result.stderr else 'N/A'}")
    except subprocess.TimeoutExpired:
        print(f"  [ERROR] acpype 超时")
        return None
    except Exception as e:
        print(f"  [ERROR] acpype 异常: {e}")
        return None
    
    # 查找生成的 itp 文件
    itp_patterns = [
        work_dir / f"{mol_name}.acpype" / f"{mol_name}_GMX.itp",
        work_dir / f"{mol_name}.acpype" / f"{mol_name}.itp",
    ]
    
    for pattern in itp_patterns:
        if pattern.exists():
            return pattern
    
    # 搜索任意 itp
    itp_files = list((work_dir).rglob("*.itp"))
    if itp_files:
        return itp_files[0]
    
    print(f"  [ERROR] 未找到 {mol_name} 的 itp 文件")
    return None


def create_simple_itp_for_ion(mol_name: str, pdb_file: Path, charge: int, output_dir: Path) -> Path:
    """为简单离子创建 itp 文件（不需要 acpype）"""
    print(f"  [simple] 为离子 {mol_name} 创建简化 itp...")
    
    # 读取 PDB 获取原子信息
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                residue = line[17:20].strip()
                atoms.append((atom_name, residue))
    
    if not atoms:
        print(f"  [ERROR] {pdb_file} 中没有原子")
        return None
    
    itp_path = output_dir / f"{mol_name}.itp"
    
    # 对于 Li+，使用标准 OPLS-AA 参数
    if mol_name.upper() == "LI":
        with open(itp_path, 'w') as f:
            f.write(f"; {mol_name} topology (Li+ ion)\n\n")
            f.write("[ moleculetype ]\n")
            f.write(f"; name  nrexcl\n")
            f.write(f"{mol_name}  3\n\n")
            f.write("[ atoms ]\n")
            f.write(";   nr  type  resnr  residue  atom  cgnr  charge   mass\n")
            f.write(f"     1  Li       1    LI      LI     1   1.000   6.94\n")
    else:
        # 对于其他离子，使用通用模板
        with open(itp_path, 'w') as f:
            f.write(f"; {mol_name} topology\n\n")
            f.write("[ moleculetype ]\n")
            f.write(f"; name  nrexcl\n")
            f.write(f"{mol_name}  3\n\n")
            f.write("[ atoms ]\n")
            f.write(";   nr  type  resnr  residue  atom  cgnr  charge   mass\n")
            
            charge_per_atom = charge / len(atoms) if atoms else 0
            for i, (atom_name, residue) in enumerate(atoms, 1):
                # 简单映射原子类型
                element = atom_name[0]
                if element == 'C':
                    atype = 'opls_135'  # 通用 C
                    mass = 12.01
                elif element == 'H':
                    atype = 'opls_140'  # 通用 H
                    mass = 1.008
                elif element == 'O':
                    atype = 'opls_180'  # 通用 O
                    mass = 16.00
                elif element == 'N':
                    atype = 'opls_237'  # 通用 N
                    mass = 14.01
                elif element == 'S':
                    atype = 'opls_222'  # 通用 S
                    mass = 32.06
                elif element == 'F':
                    atype = 'opls_965'  # 通用 F
                    mass = 19.00
                else:
                    atype = f'opls_{100+i}'
                    mass = 12.0
                
                f.write(f"     {i:3d}  {atype:8s}  1  {mol_name:4s}  {atom_name:4s}  {i:3d}  {charge_per_atom:8.4f}  {mass:.2f}\n")
    
    return itp_path


def create_topol_top(molecules: dict, itp_files: dict, output_dir: Path, box_size_nm: float = 25.0) -> Path:
    """创建主 topol.top 文件"""
    topol_path = output_dir / "topol.top"
    
    with open(topol_path, 'w') as f:
        f.write("; WSGPE Electrolyte Topology\n")
        f.write("; Generated by generate_topology.py\n\n")
        
        # 力场
        f.write("; Force field\n")
        f.write('#include "amber99sb-ildn.ff/forcefield.itp"\n\n')
        
        # 包含 itp 文件
        f.write("; Molecule topologies\n")
        for mol_id, itp_path in itp_files.items():
            if itp_path:
                # 复制 itp 到输出目录
                dst_itp = output_dir / itp_path.name
                if itp_path != dst_itp:
                    shutil.copy(itp_path, dst_itp)
                f.write(f'#include "{itp_path.name}"\n')
        f.write("\n")
        
        # 系统定义
        f.write("[ system ]\n")
        f.write("WSGPE Electrolyte\n\n")
        
        # 分子列表
        f.write("[ molecules ]\n")
        f.write("; Compound        nmols\n")
        for mol_id, info in molecules.items():
            mol_name = info['name']
            count = info['count']
            f.write(f"{mol_name:<16s} {count}\n")
    
    return topol_path


def pdb_to_gro(pdb_file: Path, gro_file: Path, box_size_nm: float = 25.0) -> bool:
    """将 PDB 转换为 GRO 格式"""
    print(f"[INFO] 转换 PDB -> GRO...")
    
    # 使用 gmx editconf
    cmd = [
        "gmx", "editconf",
        "-f", str(pdb_file),
        "-o", str(gro_file),
        "-box", str(box_size_nm), str(box_size_nm), str(box_size_nm),
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"[ERROR] gmx editconf 失败: {result.stderr}")
        return False
    
    return gro_file.exists()


def main():
    parser = argparse.ArgumentParser(description="生成 GROMACS 拓扑文件")
    parser.add_argument("-i", "--input", required=True, help="Packmol 输出目录")
    parser.add_argument("-c", "--config", required=True, help="配方 YAML 文件")
    parser.add_argument("-o", "--output", help="输出目录 (默认: 输入目录)")
    parser.add_argument("--box", type=float, default=25.0, help="盒子大小 (nm)")
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    config_path = Path(args.config)
    output_dir = Path(args.output) if args.output else input_dir
    
    # 检查输入
    gel_pdb = input_dir / "gel.pdb"
    if not gel_pdb.exists():
        print(f"[ERROR] 找不到 {gel_pdb}")
        return 1
    
    print(f"[INFO] 输入: {gel_pdb}")
    print(f"[INFO] 配方: {config_path}")
    print(f"[INFO] 输出: {output_dir}")
    
    # 解析配方
    molecules = parse_recipe(config_path)
    print(f"[INFO] 找到 {len(molecules)} 种分子")
    
    # 为每个分子生成 itp
    itp_files = {}
    base_dir = config_path.parent.parent  # 项目根目录
    
    for mol_id, info in molecules.items():
        mol_name = info['name']
        mol_file = base_dir / info['file']
        charge = info['charge']
        
        print(f"\n处理分子: {mol_name}")
        
        if not mol_file.exists():
            print(f"  [ERROR] 分子文件不存在: {mol_file}")
            continue
        
        # 对于简单离子（Li），使用简化方法
        if mol_name.upper() == "LI":
            itp = create_simple_itp_for_ion(mol_name, mol_file, charge, output_dir)
        else:
            # 尝试使用 acpype
            itp = run_acpype(mol_file, output_dir, mol_name, charge)
            
            # 如果 acpype 失败，使用简化方法
            if itp is None:
                print(f"  [FALLBACK] 使用简化 itp 生成...")
                itp = create_simple_itp_for_ion(mol_name, mol_file, charge, output_dir)
        
        if itp:
            itp_files[mol_id] = itp
            print(f"  [OK] {itp}")
    
    # 创建 topol.top
    print(f"\n[INFO] 创建 topol.top...")
    topol = create_topol_top(molecules, itp_files, output_dir, args.box)
    print(f"[OK] {topol}")
    
    # 转换 PDB 到 GRO
    gro_file = output_dir / "system.gro"
    if pdb_to_gro(gel_pdb, gro_file, args.box):
        print(f"[OK] {gro_file}")
    
    print(f"\n[SUCCESS] 拓扑生成完成!")
    print(f"  topol.top: {topol}")
    print(f"  system.gro: {gro_file}")
    print(f"\n下一步:")
    print(f"  ./scripts/run_gmx.sh -i {output_dir} -o {output_dir.parent}/gmx")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


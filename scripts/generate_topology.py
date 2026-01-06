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
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def parse_recipe(recipe_path: Path) -> dict:
    """解析配方文件，返回分子列表"""
    with open(recipe_path) as f:
        cfg = yaml.safe_load(f)
    
    components = []
    missing_counts = []
    
    if 'components' in cfg:
        components = cfg.get('components', [])
    else:
        recipe = cfg.get('recipe', {})
        components = recipe.get('components', [])
    
    if not components:
        for section in ("salt_solution", "polymer_matrix", "crosslinker", "photoinitiator"):
            entries = cfg.get(section, [])
            if not isinstance(entries, list):
                continue
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                if entry.get("type") == "salt":
                    count = entry.get("count")
                    if count is None:
                        missing_counts.append(entry.get("name", entry.get("id", "salt")))
                        continue
                    cation = entry.get("cation", {})
                    anion = entry.get("anion", {})
                    cation_count = count * int(cation.get("stoichiometry", 1))
                    anion_count = count * int(anion.get("stoichiometry", 1))
                    components.append({
                        "id": f"{entry.get('id', 'salt')}_cation",
                        "name": cation.get("name", "cation"),
                        "file": cation.get("file"),
                        "count": cation_count,
                        "charge": int(cation.get("charge", 0)),
                    })
                    components.append({
                        "id": f"{entry.get('id', 'salt')}_anion",
                        "name": anion.get("name", "anion"),
                        "file": anion.get("file"),
                        "count": anion_count,
                        "charge": int(anion.get("charge", 0)),
                    })
                    continue
                
                if entry.get("count") is None:
                    missing_counts.append(entry.get("name", entry.get("id", "component")))
                    continue
                
                components.append(entry)
    
    if missing_counts:
        missing_list = ", ".join(missing_counts)
        raise ValueError(
            f"以下组分缺少 count: {missing_list}\n"
            f"请先运行 recipe_to_counts.py 生成包含 count 的 recipe_resolved.yaml"
        )
    
    if not components:
        raise ValueError(
            "未找到可用于拓扑生成的 components。\n"
            "请使用包含 components 列表或已换算 count 的 recipe_resolved.yaml。"
        )
    
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
    """为简单离子创建 itp 文件（使用 GAFF2 原子类型）"""
    print(f"  [simple] 为离子 {mol_name} 创建简化 itp...")
    
    # 读取 PDB 获取原子信息
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) >= 78 else atom_name[:2].strip()
                atoms.append((atom_name, element))
    
    if not atoms:
        print(f"  [ERROR] {pdb_file} 中没有原子")
        return None
    
    itp_path = output_dir / f"{mol_name}.itp"
    
    gaff_map = {
        'C': 'c3',
        'H': 'hc',
        'O': 'o',
        'N': 'n3',
        'S': 's4',
        'F': 'f',
        'P': 'p5',
        'Li': 'Li',
    }
    
    gaff_params = {
        'c3': (3.39967e-01, 4.57730e-01),
        'hc': (2.64953e-01, 6.56888e-02),
        'o': (2.95992e-01, 8.78640e-01),
        'os': (3.00001e-01, 7.11280e-01),
        'n3': (3.25000e-01, 7.11280e-01),
        's4': (3.56359e-01, 1.04600e+00),
        'f': (3.11815e-01, 2.55224e-01),
        'p5': (3.74177e-01, 8.36800e-01),
        'Li': (1.2500e-01, 6.2760e-02),
    }
    
    mass_map = {
        'C': 12.01, 'H': 1.008, 'O': 16.00, 'N': 14.01,
        'S': 32.06, 'F': 19.00, 'P': 30.97, 'Li': 6.94,
    }
    
    atomtypes_set = set()
    atomtypes_lines = []
    for _, element in atoms:
        element = element if element else 'C'
        element = element.capitalize()
        atype = gaff_map.get(element, 'c3')
        if atype not in atomtypes_set:
            atomtypes_set.add(atype)
            sigma, epsilon = gaff_params.get(atype, (3.4e-01, 4.5e-01))
            mass = mass_map.get(element, 12.0)
            atomtypes_lines.append(
                f"{atype:4s}   {atype:4s}  0.00000  {mass:.3f}   A   {sigma:.5e}  {epsilon:.5e}"
            )
    
    if mol_name.upper() == "LI":
        with open(itp_path, 'w') as f:
            f.write(f"; {mol_name} topology (Li+ ion)\n\n")
            f.write("[ atomtypes ]\n")
            f.write("; name  bond_type  mass  charge  ptype  sigma  epsilon\n")
            for line in atomtypes_lines:
                f.write(f"{line}\n")
            f.write("\n")
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
            f.write("[ atomtypes ]\n")
            f.write("; name  bond_type  mass  charge  ptype  sigma  epsilon\n")
            for line in atomtypes_lines:
                f.write(f"{line}\n")
            f.write("\n")
            f.write("[ moleculetype ]\n")
            f.write(f"; name  nrexcl\n")
            f.write(f"{mol_name}  3\n\n")
            f.write("[ atoms ]\n")
            f.write(";   nr  type  resnr  residue  atom  cgnr  charge   mass\n")
            
            charge_per_atom = charge / len(atoms) if atoms else 0
            for i, (atom_name, element) in enumerate(atoms, 1):
                element = element if element else 'C'
                element = element.capitalize()
                atype = gaff_map.get(element, 'c3')
                mass = mass_map.get(element, 12.0)
                
                f.write(f"     {i:3d}  {atype:8s}  1  {mol_name:4s}  {atom_name:4s}  {i:3d}  {charge_per_atom:8.4f}  {mass:.2f}\n")
    
    return itp_path


def create_topol_top(molecules: dict, itp_files: dict, output_dir: Path, box_size_nm: float = 25.0) -> Path:
    """创建主 topol.top 文件"""
    topol_path = output_dir / "topol.top"
    
    with open(topol_path, 'w') as f:
        f.write("; WSGPE Electrolyte Topology\n")
        f.write("; Generated by generate_topology.py\n\n")
        
        # 力场设置 (GAFF2 与 acpype 输出保持一致)
        f.write("; Force field: GAFF2 (via acpype-generated itp)\n")
        f.write("[ defaults ]\n")
        f.write("; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
        f.write("1         2          yes        0.5      0.8333\n\n")
        
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
    try:
        molecules = parse_recipe(config_path)
    except ValueError as exc:
        print(f"[ERROR] {exc}")
        return 1
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

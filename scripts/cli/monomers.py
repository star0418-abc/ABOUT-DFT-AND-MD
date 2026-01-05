#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
monomers.py - 统一单体处理 CLI
==============================

合并了原 mol2_to_polymer_mol2.py 和 process_monomers_batch.py 的功能。

功能:
  - 从单体 MOL2 生成寡聚体 MOL2
  - 批量处理多个单体
  - 支持 recipe.yaml 配置

用法:
  # 批量处理（推荐）
  python -m scripts.cli.monomers --names EGDA MMA TFEMA VA

  # 单文件处理
  python -m scripts.cli.monomers --input mol2/EGDA.mol2

  # 使用 recipe 配置
  python -m scripts.cli.monomers --recipe config/recipe.yaml

  # 指定聚合度
  python -m scripts.cli.monomers --n_repeat 5
"""

import argparse
import os
import sys
import random
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path

# 添加项目根目录到 PATH
script_dir = Path(__file__).parent
project_root = script_dir.parent.parent
sys.path.insert(0, str(project_root))

# 导入公共库
from scripts.lib.recipe import load_recipe, get_oligomer_n, get_project_root
from scripts.lib.io_mol import read_mol2, write_mol2_strict, validate_mol2_input, get_atom_names_from_mol2
from scripts.lib.naming import get_polymer_output_path, get_polymer_name
from scripts.lib.validate import validate_mol2_output
from scripts.lib.logging_utils import log_success, log_error, log_info, print_header, print_summary
from scripts.lib.connectivity import validate_single_molecule, check_mol2_connectivity

# RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFHasAllMoleculeParams
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
except ImportError:
    print("[FATAL] RDKit 未安装: conda install -c conda-forge rdkit")
    sys.exit(1)


# ==============================================================================
# 常量与配置
# ==============================================================================

DEFAULT_OLIGOMER_N = 3
DEFAULT_MONOMERS = ["EGDA", "MMA", "TFEMA", "VA"]


@dataclass
class MonomerConfig:
    """单体配置"""
    name: str
    polymer_name: str
    head_atom_names: List[str]
    tail_atom_names: List[str]
    head_smarts: Optional[str]
    tail_smarts: Optional[str]
    description: str
    mw_g_mol: float = 0.0


MONOMER_CONFIGS: Dict[str, MonomerConfig] = {
    "EGDA": MonomerConfig(
        name="EGDA", polymer_name="PEGDA",
        head_atom_names=["HA", "C1"], tail_atom_names=["TB", "HB"],
        head_smarts="[CH2;$(C=C)]", tail_smarts="[CH2;$(C=C)]",
        description="乙二醇二丙烯酸酯 → PEGDA", mw_g_mol=198.22,
    ),
    "MMA": MonomerConfig(
        name="MMA", polymer_name="PMMA",
        head_atom_names=["C3", "C1"], tail_atom_names=["C2", "C4"],
        head_smarts="[CH2;$(C=C)]", tail_smarts="[C;$(C(=C)(C)C)]",
        description="甲基丙烯酸甲酯 → PMMA", mw_g_mol=100.12,
    ),
    "TFEMA": MonomerConfig(
        name="TFEMA", polymer_name="PTFEMA",
        head_atom_names=["C3", "C1"], tail_atom_names=["C2", "C4"],
        head_smarts="[CH2;$(C=C)]", tail_smarts="[C;$(C(=C)(C)C)]",
        description="三氟乙基甲基丙烯酸酯 → PTFEMA", mw_g_mol=168.12,
    ),
    "VA": MonomerConfig(
        name="VA", polymer_name="PVA",
        head_atom_names=["C6", "C4", "C1"], tail_atom_names=["C5", "C3", "C2"],
        head_smarts="[CH2;$(C=C)]", tail_smarts="[CH;$(C=C)]",
        description="乙酸乙烯酯 → PVA", mw_g_mol=86.09,
    ),
    "EO": MonomerConfig(
        name="EO", polymer_name="PEO",
        head_atom_names=["O1"], tail_atom_names=["C2"],
        head_smarts="[O;X2;R]", tail_smarts="[C;X4;R]",
        description="环氧乙烷 → PEO", mw_g_mol=44.05,
    ),
    "AM": MonomerConfig(
        name="AM", polymer_name="PAM",
        head_atom_names=["C1", "HA"], tail_atom_names=["C2", "TA"],
        head_smarts="[CH2;$(C=C)]", tail_smarts="[CH;$(C=C)]",
        description="丙烯酰胺 → PAM", mw_g_mol=71.08,
    ),
}


# ==============================================================================
# 聚合逻辑
# ==============================================================================

def find_connecting_atoms(
    mol: Chem.Mol,
    config: MonomerConfig,
    atom_names: Dict[int, str]
) -> Tuple[int, int]:
    """找到单体的 head 和 tail 连接原子"""
    head_idx = None
    tail_idx = None
    
    name_to_idx = {name: idx for idx, name in atom_names.items()}
    
    for name in config.head_atom_names:
        if name in name_to_idx:
            head_idx = name_to_idx[name]
            break
    
    for name in config.tail_atom_names:
        if name in name_to_idx:
            tail_idx = name_to_idx[name]
            break
    
    if head_idx is None and config.head_smarts:
        pattern = Chem.MolFromSmarts(config.head_smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                head_idx = matches[0][0]
    
    if tail_idx is None and config.tail_smarts:
        pattern = Chem.MolFromSmarts(config.tail_smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                for match in reversed(matches):
                    if match[0] != head_idx:
                        tail_idx = match[0]
                        break
    
    if head_idx is not None and tail_idx is None:
        pattern = Chem.MolFromSmarts("[CH2;$(C=C)]")
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if match[0] != head_idx:
                    tail_idx = match[0]
                    break
    
    if head_idx is None or tail_idx is None:
        raise RuntimeError(f"无法识别连接位点: Head={head_idx}, Tail={tail_idx}")
    
    return head_idx, tail_idx


def find_sacrificial_hydrogen(mol: Chem.Mol, atom_idx: int) -> Optional[int]:
    """找到连接原子相邻的牺牲氢"""
    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            return neighbor.GetIdx()
    return None


def prepare_monomer(mol: Chem.Mol, head_idx: int, tail_idx: int) -> Chem.Mol:
    """准备单体用于聚合"""
    mol = Chem.RWMol(Chem.Mol(mol))
    bond = mol.GetBondBetweenAtoms(head_idx, tail_idx)
    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
        bond.SetBondType(Chem.BondType.SINGLE)
        mol.GetAtomWithIdx(head_idx).SetNumExplicitHs(
            mol.GetAtomWithIdx(head_idx).GetNumExplicitHs() + 1)
        mol.GetAtomWithIdx(tail_idx).SetNumExplicitHs(
            mol.GetAtomWithIdx(tail_idx).GetNumExplicitHs() + 1)
    return mol.GetMol()


def polymerize(
    monomer: Chem.Mol,
    head_idx: int,
    tail_idx: int,
    dp: int,
    seed: int = 2025
) -> Chem.Mol:
    """生成线性聚合物"""
    if dp < 1:
        raise ValueError("聚合度必须 >= 1")
    if dp == 1:
        return Chem.Mol(monomer)
    
    random.seed(seed)
    prep_monomer = prepare_monomer(monomer, head_idx, tail_idx)
    prep_monomer = Chem.AddHs(prep_monomer, addCoords=True)
    
    polymer = Chem.RWMol(Chem.Mol(prep_monomer))
    current_tail = tail_idx
    
    for i in range(1, dp):
        offset = polymer.GetNumAtoms()
        combined = Chem.CombineMols(polymer.GetMol(), Chem.Mol(prep_monomer))
        polymer = Chem.RWMol(combined)
        
        new_head = head_idx + offset
        new_tail = tail_idx + offset
        
        tail_h = find_sacrificial_hydrogen(polymer, current_tail)
        head_h = find_sacrificial_hydrogen(polymer, new_head)
        
        hydrogens_to_remove = sorted(
            [h for h in [tail_h, head_h] if h is not None],
            reverse=True
        )
        
        tail_adjust = sum(1 for h in hydrogens_to_remove if h < current_tail)
        head_adjust = sum(1 for h in hydrogens_to_remove if h < new_head)
        new_tail_adjust = sum(1 for h in hydrogens_to_remove if h < new_tail)
        
        for h in hydrogens_to_remove:
            polymer.RemoveAtom(h)
        
        current_tail -= tail_adjust
        new_head -= head_adjust
        new_tail -= new_tail_adjust
        
        polymer.AddBond(current_tail, new_head, Chem.BondType.SINGLE)
        current_tail = new_tail
    
    result = polymer.GetMol()
    try:
        Chem.SanitizeMol(result)
    except:
        pass
    
    try:
        result.RemoveAllConformers()
        status = AllChem.EmbedMolecule(result, maxAttempts=50, randomSeed=seed)
        if status < 0:
            AllChem.EmbedMolecule(result, useRandomCoords=True, randomSeed=seed)
    except:
        result.RemoveAllConformers()
        AllChem.EmbedMolecule(result, useRandomCoords=True, randomSeed=seed, maxAttempts=100)
    
    return result


def process_single_monomer(
    input_path: str,
    oligomer_n: int,
    optimize: bool = True,
    seed: int = 2025
) -> str:
    """
    处理单个单体文件
    
    Returns:
        输出文件路径
    """
    validate_mol2_input(input_path)
    
    input_path = os.path.abspath(input_path)
    output_path = get_polymer_output_path(input_path)
    input_name = os.path.splitext(os.path.basename(input_path))[0].upper()
    
    config = MONOMER_CONFIGS.get(input_name, MonomerConfig(
        name=input_name, polymer_name=get_polymer_name(input_name),
        head_atom_names=["HA", "C1", "C3"], tail_atom_names=["TA", "TB", "C2"],
        head_smarts="[CH2;$(C=C)]", tail_smarts="[CH;$(C=C)]",
        description=f"通用聚合: {input_name}",
    ))
    
    print_header(f"处理: {config.name} → {config.polymer_name} (n={oligomer_n})")
    print(f"  输入路径: {input_path}")
    print(f"  输出路径: {output_path}")
    print(f"  聚合度:   {oligomer_n}")
    
    monomer = read_mol2(input_path)
    monomer_atoms = monomer.GetNumAtoms()
    monomer_bonds = monomer.GetNumBonds()
    print(f"  输入原子: {monomer_atoms}")
    print(f"  输入键数: {monomer_bonds}")
    
    atom_names = get_atom_names_from_mol2(input_path)
    head_idx, tail_idx = find_connecting_atoms(monomer, config, atom_names)
    print(f"  连接位点: Head={head_idx + 1}, Tail={tail_idx + 1}")
    
    polymer = polymerize(monomer, head_idx, tail_idx, oligomer_n, seed)
    polymer_atoms = polymer.GetNumAtoms()
    polymer_bonds = polymer.GetNumBonds()
    print(f"  输出原子: {polymer_atoms}")
    print(f"  输出键数: {polymer_bonds}")
    
    if optimize:
        print("  UFF 优化...")
        try:
            if UFFHasAllMoleculeParams(polymer):
                UFFOptimizeMolecule(polymer, maxIters=500)
        except:
            pass
    
    write_mol2_strict(polymer, output_path, config.polymer_name, config.polymer_name[:3])
    
    size = validate_mol2_output(output_path)
    
    # 连通性检查（强制）
    validate_single_molecule(output_path)
    
    # 验证聚合物片段数
    frags = Chem.GetMolFrags(polymer)
    if len(frags) != 1:
        raise RuntimeError(f"聚合物不是单一连通分子！发现 {len(frags)} 个碎片")
    
    print("")
    log_success(f"写文件成功: {output_path} ({size} bytes)")
    
    return output_path


def process_batch(
    mol2_dir: str,
    monomers: List[str],
    oligomer_n: int,
    optimize: bool = True,
    seed: int = 2025
) -> Dict[str, List]:
    """批量处理"""
    results = {"success": [], "failed": []}
    
    print_header("单体批处理: MOL2 → 寡聚体 MOL2")
    print(f"  目录:   {mol2_dir}")
    print(f"  单体:   {', '.join(monomers)}")
    print(f"  聚合度: {oligomer_n}")
    
    for monomer in monomers:
        input_path = os.path.join(mol2_dir, f"{monomer}.mol2")
        
        if not os.path.isfile(input_path):
            log_error(f"单体文件不存在: {input_path}")
            results["failed"].append(monomer)
            continue
        
        try:
            output_path = process_single_monomer(input_path, oligomer_n, optimize, seed)
            results["success"].append((monomer, output_path))
        except Exception as e:
            log_error(f"处理 {monomer} 失败: {e}")
            results["failed"].append(monomer)
    
    print_summary("批处理结果", {
        "成功": len(results["success"]),
        "失败": len(results["failed"]),
    })
    
    for name, path in results["success"]:
        log_success(f"{name} → {os.path.basename(path)}")
    for name in results["failed"]:
        log_error(name)
    
    return results


# ==============================================================================
# CLI
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="单体 MOL2 → 寡聚体 MOL2 (默认三聚体)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
示例:
  # 批量处理四个单体
  python -m scripts.cli.monomers --names EGDA MMA TFEMA VA

  # 单文件处理
  python -m scripts.cli.monomers --input mol2/EGDA.mol2

  # 使用 recipe 配置
  python -m scripts.cli.monomers --recipe config/recipe.yaml

默认单体: {', '.join(DEFAULT_MONOMERS)}
默认聚合度: {DEFAULT_OLIGOMER_N}
"""
    )
    
    parser.add_argument("--input", "-i", help="单个输入 MOL2 文件")
    parser.add_argument("--dir", default="mol2", help="MOL2 目录 (批处理)")
    parser.add_argument("--names", nargs="+", default=None, help="单体列表 (批处理)")
    parser.add_argument("--n_repeat", type=int, default=None, help="聚合度")
    parser.add_argument("--recipe", "-r", default=None, help="recipe.yaml 路径")
    parser.add_argument("--no_opt", action="store_true", help="禁用 UFF 优化")
    parser.add_argument("--seed", type=int, default=2025, help="随机种子")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    
    # 确定 oligomer_n
    config = None
    if args.recipe:
        try:
            config = load_recipe(args.recipe)
        except Exception as e:
            log_error(f"加载 recipe 失败: {e}")
    
    oligomer_n = get_oligomer_n(config, args.n_repeat, DEFAULT_OLIGOMER_N)
    source = "CLI" if args.n_repeat else ("recipe" if config else "默认")
    log_info(f"聚合度: {oligomer_n} (来源: {source})")
    
    if args.input:
        # 单文件模式
        try:
            process_single_monomer(args.input, oligomer_n, not args.no_opt, args.seed)
            return 0
        except Exception as e:
            log_error(str(e))
            return 1
    else:
        # 批处理模式
        monomers = [m.upper() for m in (args.names or DEFAULT_MONOMERS)]
        mol2_dir = os.path.join(get_project_root(), args.dir)
        
        if not os.path.isdir(mol2_dir):
            log_error(f"目录不存在: {mol2_dir}")
            return 1
        
        results = process_batch(mol2_dir, monomers, oligomer_n, not args.no_opt, args.seed)
        return 1 if results["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())


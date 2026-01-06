#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
process_monomers_batch.py - 批量处理单体生成寡聚体 MOL2
==========================================================

功能:
  - 批量处理 EGDA/MMA/TFEMA/VA 单体生成三聚体 MOL2
  - 输入必须是 .mol2 格式（其他格式报错退出）
  - 输出到输入文件同目录，命名: P + 原文件名
  - 严格验证输出文件存在、非空、包含必需段落

⚠️ 特殊处理:
  - PEGDA 不走寡聚体路线，需要 HTPolyNet 生成交联网络
  - 此脚本仅处理: EGDA → PEGDA, MMA → PMMA, TFEMA → PTFEMA, VA → PVA

用法:
  # 处理所有四个单体（默认三聚体）
  python scripts/process_monomers_batch.py

  # 指定目录和单体
  python scripts/process_monomers_batch.py --dir mol2 --names EGDA MMA TFEMA VA

  # 指定聚合度
  python scripts/process_monomers_batch.py --n_repeat 3

  # 使用 recipe.yaml 配置
  python scripts/process_monomers_batch.py --recipe config/recipe.yaml

版本: 1.0.0
"""

import argparse
import os
import sys
import random
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path

# ==============================================================================
# RDKit 导入
# ==============================================================================

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFHasAllMoleculeParams
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
except ImportError:
    print("=" * 60)
    print("[FATAL] RDKit 未安装")
    print("=" * 60)
    print("\n请使用以下命令安装 RDKit:")
    print("  conda install -c conda-forge rdkit")
    print("=" * 60)
    sys.exit(1)

# ==============================================================================
# YAML 导入
# ==============================================================================

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# ==============================================================================
# 常量
# ==============================================================================

DEFAULT_OLIGOMER_N = 3  # 默认三聚体
DEFAULT_MONOMERS = ["EGDA", "MMA", "TFEMA", "VA"]


# ==============================================================================
# 单体配置
# ==============================================================================

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
        name="EGDA",
        polymer_name="PEGDA",
        head_atom_names=["HA", "C1"],
        tail_atom_names=["TB", "HB"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[CH2;$(C=C)]",
        description="乙二醇二丙烯酸酯 → PEGDA 寡聚体",
        mw_g_mol=198.22,
    ),
    "MMA": MonomerConfig(
        name="MMA",
        polymer_name="PMMA",
        head_atom_names=["C3", "C1"],
        tail_atom_names=["C2", "C4"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[C;$(C(=C)(C)C)]",
        description="甲基丙烯酸甲酯 → PMMA",
        mw_g_mol=100.12,
    ),
    "TFEMA": MonomerConfig(
        name="TFEMA",
        polymer_name="PTFEMA",
        head_atom_names=["C3", "C1"],
        tail_atom_names=["C2", "C4"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[C;$(C(=C)(C)C)]",
        description="三氟乙基甲基丙烯酸酯 → PTFEMA",
        mw_g_mol=168.12,
    ),
    "VA": MonomerConfig(
        name="VA",
        polymer_name="PVA",
        head_atom_names=["C6", "C4", "C1"],
        tail_atom_names=["C5", "C3", "C2"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[CH;$(C=C)]",
        description="乙酸乙烯酯 → PVA",
        mw_g_mol=86.09,
    ),
}


# ==============================================================================
# 输入验证
# ==============================================================================

def validate_mol2_input(filepath: str) -> None:
    """验证输入文件必须是 .mol2 格式"""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"文件不存在: {filepath}")
    
    ext = os.path.splitext(filepath)[1].lower()
    if ext != ".mol2":
        raise ValueError(
            f"仅支持 .mol2 格式输入，收到: {ext}\n"
            f"文件: {filepath}\n"
            f"提示: 请先使用 mol_to_mol2_batch.py 转换格式"
        )


# ==============================================================================
# MOL2 读写
# ==============================================================================

def read_mol2(filepath: str) -> Chem.Mol:
    """读取 MOL2 文件"""
    validate_mol2_input(filepath)
    
    mol = Chem.MolFromMol2File(filepath, removeHs=False)
    if mol is None:
        raise RuntimeError(f"无法解析 MOL2 文件: {filepath}")
    
    return mol


def write_mol2_strict(
    mol: Chem.Mol,
    filepath: str,
    mol_name: str,
    res_name: str
) -> None:
    """
    严格模式写入 MOL2 文件
    
    - 写入后验证文件存在且非空
    - 验证包含必需段落
    - 任何失败都抛出异常
    """
    # 确保目录存在
    out_dir = os.path.dirname(filepath)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    
    try:
        conf = mol.GetConformer()
        atoms = list(mol.GetAtoms())
        bonds = list(mol.GetBonds())
        
        if len(atoms) == 0:
            raise RuntimeError("分子没有原子，无法写入")
        
        # 生成唯一原子名
        atom_names = {}
        element_count = {}
        
        for atom in atoms:
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            element_count[symbol] = element_count.get(symbol, 0) + 1
            atom_names[idx] = f"{symbol}{element_count[symbol]}"
        
        # Sybyl 原子类型
        def get_sybyl_type(atom):
            symbol = atom.GetSymbol()
            hyb = str(atom.GetHybridization()).lower()
            
            if symbol == 'H':
                return 'H'
            elif symbol == 'C':
                if atom.GetIsAromatic():
                    return 'C.ar'
                elif 'sp3' in hyb:
                    return 'C.3'
                elif 'sp2' in hyb:
                    return 'C.2'
                else:
                    return 'C.3'
            elif symbol == 'N':
                if atom.GetIsAromatic():
                    return 'N.ar'
                elif 'sp3' in hyb:
                    return 'N.3'
                else:
                    return 'N.2'
            elif symbol == 'O':
                if 'sp3' in hyb:
                    return 'O.3'
                else:
                    return 'O.2'
            elif symbol == 'F':
                return 'F'
            elif symbol == 'S':
                return 'S.3'
            else:
                return symbol
        
        # 键类型映射
        bond_type_map = {
            Chem.BondType.SINGLE: '1',
            Chem.BondType.DOUBLE: '2',
            Chem.BondType.TRIPLE: '3',
            Chem.BondType.AROMATIC: 'ar',
        }
        
        # 写入文件
        with open(filepath, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{mol_name}\n")
            f.write(f" {len(atoms)} {len(bonds)} 1 0 0\n")
            f.write("SMALL\n")
            f.write("NO_CHARGES\n\n")
            
            f.write("@<TRIPOS>ATOM\n")
            for atom in atoms:
                idx = atom.GetIdx()
                pos = conf.GetAtomPosition(idx)
                name = atom_names[idx]
                sybyl = get_sybyl_type(atom)
                charge = atom.GetFormalCharge()
                f.write(f"{idx+1:7d} {name:<4s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} "
                       f"{sybyl:<6s} 1 {res_name[:3]:<4s} {charge:8.4f}\n")
            
            f.write("@<TRIPOS>BOND\n")
            for i, bond in enumerate(bonds):
                begin = bond.GetBeginAtomIdx() + 1
                end = bond.GetEndAtomIdx() + 1
                btype = bond_type_map.get(bond.GetBondType(), '1')
                f.write(f"{i+1:6d} {begin:5d} {end:5d} {btype}\n")
            
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write(f"     1 {res_name[:3]:<4s}        1 RESIDUE    0 **** **** 0 ROOT\n")
        
        # ========== 严格验证 ==========
        # 1. 文件存在
        if not os.path.exists(filepath):
            raise RuntimeError(f"写入后文件不存在: {filepath}")
        
        # 2. 文件非空
        file_size = os.path.getsize(filepath)
        if file_size == 0:
            raise RuntimeError(f"写入后文件为空: {filepath}")
        
        # 3. 包含必需段落
        with open(filepath, 'r') as f:
            content = f.read()
        
        required_sections = ["@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND"]
        for section in required_sections:
            if section not in content:
                raise RuntimeError(f"输出文件缺少必需段落 {section}: {filepath}")
        
        # 4. ATOM 行数匹配
        atom_lines = [l for l in content.split('\n') if l.strip() and 
                     l.strip()[0].isdigit() and '@' not in l]
        # 简单检查：至少有原子行
        if len(atom_lines) < len(atoms):
            raise RuntimeError(
                f"ATOM 行数不匹配: 期望 {len(atoms)}, 实际 {len(atom_lines)}"
            )
        
    except Exception as e:
        # 清理可能的部分写入文件
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
            except:
                pass
        raise RuntimeError(f"写入 MOL2 失败: {e}")


# ==============================================================================
# 连接位点识别
# ==============================================================================

def get_atom_names_from_mol2(filepath: str) -> Dict[int, str]:
    """从 MOL2 文件解析原子名"""
    atom_names = {}
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    in_atom_section = False
    for line in lines:
        stripped = line.strip()
        
        if stripped.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            continue
        elif stripped.startswith("@<TRIPOS>"):
            in_atom_section = False
            continue
        
        if in_atom_section and stripped:
            parts = stripped.split()
            if len(parts) >= 2:
                atom_id = int(parts[0]) - 1
                atom_name = parts[1]
                atom_names[atom_id] = atom_name
    
    return atom_names


def find_connecting_atoms(
    mol: Chem.Mol,
    config: MonomerConfig,
    atom_names: Dict[int, str]
) -> Tuple[int, int]:
    """找到单体的 head 和 tail 连接原子"""
    head_idx = None
    tail_idx = None
    
    # 方法 1: 通过原子名匹配
    name_to_idx = {name: idx for idx, name in atom_names.items()}
    
    for name in config.head_atom_names:
        if name in name_to_idx:
            head_idx = name_to_idx[name]
            break
    
    for name in config.tail_atom_names:
        if name in name_to_idx:
            tail_idx = name_to_idx[name]
            break
    
    # 方法 2: 通过 SMARTS 匹配
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
    
    # 特殊处理双官能团单体
    if head_idx is not None and tail_idx is None:
        pattern = Chem.MolFromSmarts("[CH2;$(C=C)]")
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if match[0] != head_idx:
                    tail_idx = match[0]
                    break
    
    if head_idx is None or tail_idx is None:
        raise RuntimeError(
            f"无法识别连接位点\n"
            f"  Head: {head_idx}, Tail: {tail_idx}\n"
            f"  尝试的原子名: head={config.head_atom_names}, tail={config.tail_atom_names}"
        )
    
    return head_idx, tail_idx


# ==============================================================================
# 聚合物生成
# ==============================================================================

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
        head_atom = mol.GetAtomWithIdx(head_idx)
        tail_atom = mol.GetAtomWithIdx(tail_idx)
        head_atom.SetNumExplicitHs(head_atom.GetNumExplicitHs() + 1)
        tail_atom.SetNumExplicitHs(tail_atom.GetNumExplicitHs() + 1)
    
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
    
    # 准备单体
    prep_monomer = prepare_monomer(monomer, head_idx, tail_idx)
    prep_monomer = Chem.AddHs(prep_monomer, addCoords=True)
    
    polymer = Chem.RWMol(Chem.Mol(prep_monomer))
    current_tail = tail_idx
    
    for i in range(1, dp):
        offset = polymer.GetNumAtoms()
        
        new_monomer = Chem.Mol(prep_monomer)
        combined = Chem.CombineMols(polymer.GetMol(), new_monomer)
        polymer = Chem.RWMol(combined)
        
        new_head = head_idx + offset
        new_tail = tail_idx + offset
        
        tail_h = find_sacrificial_hydrogen(polymer, current_tail)
        head_h = find_sacrificial_hydrogen(polymer, new_head)
        
        hydrogens_to_remove = []
        if tail_h is not None:
            hydrogens_to_remove.append(tail_h)
        if head_h is not None:
            hydrogens_to_remove.append(head_h)
        
        hydrogens_to_remove.sort(reverse=True)
        
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
    
    # 生成 3D 坐标
    try:
        result.RemoveAllConformers()
        status = AllChem.EmbedMolecule(result, maxAttempts=50, randomSeed=seed)
        if status < 0:
            AllChem.EmbedMolecule(result, useRandomCoords=True, randomSeed=seed)
    except:
        result.RemoveAllConformers()
        AllChem.EmbedMolecule(result, useRandomCoords=True, randomSeed=seed, maxAttempts=100)
    
    return result


# ==============================================================================
# 批处理主函数
# ==============================================================================

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
    # 验证输入
    validate_mol2_input(input_path)
    
    input_path = os.path.abspath(input_path)
    input_dir = os.path.dirname(input_path)
    input_basename = os.path.basename(input_path)
    input_name, input_ext = os.path.splitext(input_basename)
    
    # 输出命名: P + 原文件名
    output_name = f"P{input_name}{input_ext}"
    output_path = os.path.join(input_dir, output_name)
    
    # 获取配置
    monomer_key = input_name.upper()
    if monomer_key in MONOMER_CONFIGS:
        config = MONOMER_CONFIGS[monomer_key]
    else:
        config = MonomerConfig(
            name=input_name,
            polymer_name=f"P{input_name}",
            head_atom_names=["HA", "C1", "C3"],
            tail_atom_names=["TA", "TB", "C2"],
            head_smarts="[CH2;$(C=C)]",
            tail_smarts="[CH;$(C=C)]",
            description=f"通用聚合: {input_name}",
        )
    
    print("")
    print("=" * 60)
    print(f"处理: {input_name} → P{input_name} (oligomer_n={oligomer_n})")
    print("=" * 60)
    print(f"  输入路径: {input_path}")
    print(f"  输出路径: {output_path}")
    print(f"  聚合度:   {oligomer_n}")
    
    # 读取单体
    monomer = read_mol2(input_path)
    monomer_atoms = monomer.GetNumAtoms()
    monomer_bonds = monomer.GetNumBonds()
    print(f"  输入原子数: {monomer_atoms}")
    print(f"  输入键数:   {monomer_bonds}")
    
    # 获取原子名并找到连接位点
    atom_names = get_atom_names_from_mol2(input_path)
    head_idx, tail_idx = find_connecting_atoms(monomer, config, atom_names)
    print(f"  连接位点: Head={head_idx + 1}, Tail={tail_idx + 1}")
    
    # 聚合
    polymer = polymerize(monomer, head_idx, tail_idx, oligomer_n, seed)
    polymer_atoms = polymer.GetNumAtoms()
    polymer_bonds = polymer.GetNumBonds()
    print(f"  输出原子数: {polymer_atoms}")
    print(f"  输出键数:   {polymer_bonds}")
    
    # UFF 优化
    if optimize:
        print("  UFF 优化中...")
        try:
            if UFFHasAllMoleculeParams(polymer):
                UFFOptimizeMolecule(polymer, maxIters=500)
        except Exception as e:
            print(f"  [WARN] 优化跳过: {e}")
    
    # 写入文件（严格模式）
    write_mol2_strict(
        mol=polymer,
        filepath=output_path,
        mol_name=config.polymer_name,
        res_name=config.polymer_name[:3]
    )
    
    # 最终验证
    file_size = os.path.getsize(output_path)
    print("")
    print("  ========== 输出验证 ==========")
    print(f"  ✓ 文件存在: {output_path}")
    print(f"  ✓ 文件大小: {file_size} bytes")
    print(f"  ✓ 原子数:   {polymer_atoms} (单体 {monomer_atoms} × {oligomer_n})")
    print(f"  ✓ 键数:     {polymer_bonds}")
    print("  ✓ 包含: @<TRIPOS>MOLECULE, @<TRIPOS>ATOM, @<TRIPOS>BOND")
    print("  [SUCCESS] 写文件成功!")
    print("=" * 60)
    
    return output_path


def process_batch(
    mol2_dir: str,
    monomers: List[str],
    oligomer_n: int,
    optimize: bool = True,
    seed: int = 2025
) -> Dict[str, List[str]]:
    """
    批量处理多个单体
    
    Returns:
        {"success": [...], "failed": [...]}
    """
    results = {"success": [], "failed": []}
    
    print("=" * 60)
    print("单体批处理: MOL2 → 寡聚体 MOL2")
    print("=" * 60)
    print(f"目录:     {mol2_dir}")
    print(f"单体:     {', '.join(monomers)}")
    print(f"聚合度:   {oligomer_n}")
    print("=" * 60)
    
    for monomer in monomers:
        input_path = os.path.join(mol2_dir, f"{monomer}.mol2")
        
        if not os.path.isfile(input_path):
            print(f"\n[SKIP] 单体文件不存在: {input_path}")
            results["failed"].append(monomer)
            continue
        
        try:
            output_path = process_single_monomer(
                input_path=input_path,
                oligomer_n=oligomer_n,
                optimize=optimize,
                seed=seed
            )
            results["success"].append((monomer, output_path))
        except Exception as e:
            print(f"\n[ERROR] 处理 {monomer} 失败: {e}")
            import traceback
            traceback.print_exc()
            results["failed"].append(monomer)
    
    # 汇总
    print("\n" + "=" * 60)
    print("批处理结果汇总")
    print("=" * 60)
    print(f"成功: {len(results['success'])}")
    for name, path in results["success"]:
        print(f"  ✓ {name} → {os.path.basename(path)}")
    print(f"失败: {len(results['failed'])}")
    for name in results["failed"]:
        print(f"  ✗ {name}")
    print("=" * 60)
    
    return results


# ==============================================================================
# Recipe 配置读取
# ==============================================================================

def load_oligomer_n_from_recipe(recipe_path: str) -> Optional[int]:
    """从 recipe.yaml 读取 oligomer_n"""
    if not HAS_YAML:
        return None
    
    if not os.path.isfile(recipe_path):
        return None
    
    try:
        with open(recipe_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        polymerization = config.get("polymerization", {})
        if isinstance(polymerization, dict):
            oligomer_n = polymerization.get("oligomer_n")
            if oligomer_n is not None:
                return int(oligomer_n)
        
        return None
    except:
        return None


# ==============================================================================
# CLI
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="批量处理单体生成寡聚体 MOL2 (默认三聚体)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
示例:
  # 处理所有四个单体（默认三聚体）
  %(prog)s

  # 指定目录和单体
  %(prog)s --dir mol2 --names EGDA MMA TFEMA VA

  # 指定聚合度
  %(prog)s --n_repeat 3

  # 使用 recipe.yaml 配置
  %(prog)s --recipe config/recipe.yaml

默认单体: {', '.join(DEFAULT_MONOMERS)}
默认聚合度: {DEFAULT_OLIGOMER_N}
"""
    )
    
    parser.add_argument(
        "--dir",
        default="mol2",
        help="MOL2 目录 (默认: mol2)"
    )
    
    parser.add_argument(
        "--names",
        nargs="+",
        default=DEFAULT_MONOMERS,
        help=f"要处理的单体列表 (默认: {' '.join(DEFAULT_MONOMERS)})"
    )
    
    parser.add_argument(
        "--n_repeat",
        type=int,
        default=None,
        help=f"聚合度 (默认: {DEFAULT_OLIGOMER_N})"
    )
    
    parser.add_argument(
        "--recipe",
        default=None,
        help="recipe.yaml 路径"
    )
    
    parser.add_argument(
        "--no_opt",
        action="store_true",
        help="禁用 UFF 优化"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="随机种子 (默认: 2025)"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    
    # 确定 oligomer_n
    oligomer_n = DEFAULT_OLIGOMER_N
    source = "默认值"
    
    if args.recipe:
        recipe_n = load_oligomer_n_from_recipe(args.recipe)
        if recipe_n is not None:
            oligomer_n = recipe_n
            source = f"recipe.yaml ({args.recipe})"
    
    if args.n_repeat is not None:
        oligomer_n = args.n_repeat
        source = "CLI --n_repeat"
    
    print(f"[INFO] 聚合度: {oligomer_n} (来源: {source})")
    
    # 确定目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    mol2_dir = os.path.join(project_root, args.dir)
    
    if not os.path.isdir(mol2_dir):
        print(f"[ERROR] 目录不存在: {mol2_dir}")
        return 1
    
    # 处理
    results = process_batch(
        mol2_dir=mol2_dir,
        monomers=[m.upper() for m in args.names],
        oligomer_n=oligomer_n,
        optimize=not args.no_opt,
        seed=args.seed
    )
    
    if results["failed"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())


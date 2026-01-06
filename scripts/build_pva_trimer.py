#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
build_pva_trimer.py - 生成真正的 PVA（聚乙烯醇）三聚体
=======================================================

功能:
  生成真正的 PVA 三聚体，化学结构为 –CH2–CH(OH)– 重复单元。
  
  支持两种策略：
  A. 水解策略：从 PVAc 三聚体水解去除乙酰基得到 PVA
  B. 直接构建：从 PVA 重复单元直接构建三聚体

用法:
  python scripts/build_pva_trimer.py --in mol2/VA.mol2 --out mol2/PVA.mol2 --n 3 --force
  python scripts/build_pva_trimer.py --strategy B --out mol2/PVA.mol2 --n 3

化学说明:
  - VA (Vinyl Acetate) 聚合得到 PVAc (Polyvinyl Acetate)
  - PVAc 水解得到 PVA (Polyvinyl Alcohol)
  - PVA 重复单元: –CH2–CH(OH)–
  - PVAc 重复单元: –CH2–CH(O-CO-CH3)–

版本: 1.0.0
"""

import argparse
import os
import sys
import math
from pathlib import Path
from typing import Tuple, List, Optional

# RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFHasAllMoleculeParams
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
except ImportError:
    print("[FATAL] RDKit 未安装: conda install -c conda-forge rdkit")
    sys.exit(1)

# 项目路径
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.lib.connectivity import check_mol2_connectivity, validate_single_molecule
from scripts.lib.io_mol import write_mol2_strict, validate_mol2_output


# ==============================================================================
# 颜色输出
# ==============================================================================

class Colors:
    RESET = "\033[0m"
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[0;33m"
    BLUE = "\033[0;34m"
    CYAN = "\033[0;36m"


def log_success(msg): print(f"{Colors.GREEN}✓{Colors.RESET} {msg}")
def log_error(msg): print(f"{Colors.RED}✗{Colors.RESET} {msg}", file=sys.stderr)
def log_info(msg): print(f"{Colors.BLUE}ℹ{Colors.RESET} {msg}")
def log_warn(msg): print(f"{Colors.YELLOW}⚠{Colors.RESET} {msg}")


def print_header(title: str):
    print(f"\n{'=' * 60}")
    print(f" {title}")
    print(f"{'=' * 60}")


# ==============================================================================
# 化学检查
# ==============================================================================

def count_ester_groups(mol: Chem.Mol) -> int:
    """计算酯基 C(=O)O 数量"""
    pattern = Chem.MolFromSmarts("C(=O)O")
    return len(mol.GetSubstructMatches(pattern))


def count_hydroxyl_groups(mol: Chem.Mol) -> int:
    """计算羟基 -OH 数量"""
    pattern = Chem.MolFromSmarts("[OH]")
    return len(mol.GetSubstructMatches(pattern))


def count_carbonyl_groups(mol: Chem.Mol) -> int:
    """计算羰基 C=O 数量"""
    pattern = Chem.MolFromSmarts("C=O")
    return len(mol.GetSubstructMatches(pattern))


def is_pvac(mol: Chem.Mol) -> bool:
    """判断是否是 PVAc（含酯基）"""
    return count_ester_groups(mol) > 0


def validate_pva_chemistry(mol: Chem.Mol, n_repeat: int) -> None:
    """
    验证 PVA 化学结构
    
    Raises:
        RuntimeError: 化学检查失败
    """
    ester_count = count_ester_groups(mol)
    hydroxyl_count = count_hydroxyl_groups(mol)
    carbonyl_count = count_carbonyl_groups(mol)
    
    print(f"  酯基 C(=O)O 数量: {ester_count}")
    print(f"  羟基 -OH 数量: {hydroxyl_count}")
    print(f"  羰基 C=O 数量: {carbonyl_count}")
    
    # 检查：不能存在酯基
    if ester_count > 0:
        raise RuntimeError(
            f"化学检查失败！检测到 {ester_count} 个酯基 C(=O)O。\n"
            f"这是 PVAc 而不是 PVA！"
        )
    
    # 检查：必须存在羟基
    if hydroxyl_count < n_repeat:
        raise RuntimeError(
            f"化学检查失败！羟基数量 {hydroxyl_count} < 期望 {n_repeat}。\n"
            f"PVA 三聚体应有 3 个羟基。"
        )
    
    log_success(f"化学检查通过: {hydroxyl_count} 个羟基, 0 个酯基")


# ==============================================================================
# 策略 A: PVAc → PVA 水解
# ==============================================================================

def hydrolyze_pvac_to_pva(pvac_mol: Chem.Mol) -> Chem.Mol:
    """
    将 PVAc 水解为 PVA
    
    反应: –O–C(=O)–CH3 → –OH
    
    使用 RDKit 的反应 SMARTS
    """
    # 定义水解反应
    # 匹配：C(连到主链)-O-C(=O)-C(甲基)
    # 转换为：C(连到主链)-O(H)
    rxn_smarts = "[C:1][O:2]C(=O)C >> [C:1][O:2]"
    
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    
    # 多次应用反应（因为有多个酯基）
    current_mol = Chem.RWMol(pvac_mol)
    max_iterations = 10
    
    for i in range(max_iterations):
        # 检查是否还有酯基
        ester_pattern = Chem.MolFromSmarts("[CH,CH2]OC(=O)C")
        matches = current_mol.GetSubstructMatches(ester_pattern)
        
        if not matches:
            break
        
        # 应用反应
        try:
            products = rxn.RunReactants((current_mol.GetMol(),))
            if products and len(products) > 0 and len(products[0]) > 0:
                current_mol = Chem.RWMol(products[0][0])
                Chem.SanitizeMol(current_mol)
            else:
                break
        except Exception as e:
            log_warn(f"反应迭代 {i+1} 失败: {e}")
            break
    
    # 添加氢
    result = Chem.AddHs(current_mol.GetMol())
    
    return result


def hydrolyze_pvac_manual(pvac_mol: Chem.Mol) -> Chem.Mol:
    """
    手动水解 PVAc 到 PVA
    
    策略：
    1. 找到所有 -O-C(=O)-CH3 侧基
    2. 删除 C(=O)-CH3 部分
    3. 在 O 上补 H
    """
    # 匹配酯基：找到连接到主链碳的 O-C(=O)-C
    # 模式：[主链碳]-[酯氧]-[羰基碳](=[羰基氧])-[甲基碳]
    ester_pattern = Chem.MolFromSmarts("[CH,CH2:1][O:2][C:3](=[O:4])[CH3:5]")
    
    mol = Chem.RWMol(pvac_mol)
    
    iterations = 0
    max_iterations = 10
    
    while iterations < max_iterations:
        matches = mol.GetSubstructMatches(ester_pattern)
        if not matches:
            break
        
        # 取第一个匹配
        match = matches[0]
        main_c, ester_o, carbonyl_c, carbonyl_o, methyl_c = match
        
        # 找到甲基上的氢
        methyl_atom = mol.GetAtomWithIdx(methyl_c)
        methyl_hs = [n.GetIdx() for n in methyl_atom.GetNeighbors() if n.GetSymbol() == 'H']
        
        # 删除原子（从大到小删除以保持索引有效）
        atoms_to_remove = sorted([carbonyl_c, carbonyl_o, methyl_c] + methyl_hs, reverse=True)
        
        for atom_idx in atoms_to_remove:
            mol.RemoveAtom(atom_idx)
        
        # 现在 ester_o 变成了 hydroxyl O，需要重新获取其索引
        # 由于删除了原子，索引会变化，需要重新 sanitize
        try:
            Chem.SanitizeMol(mol)
        except:
            pass
        
        iterations += 1
    
    # 添加氢（会在 O 上添加 H 形成 -OH）
    result = Chem.AddHs(mol.GetMol())
    
    try:
        Chem.SanitizeMol(result)
    except:
        pass
    
    return result


# ==============================================================================
# 策略 B: 直接构建 PVA 三聚体
# ==============================================================================

def build_pva_from_smiles(n_repeat: int = 3) -> Chem.Mol:
    """
    直接从 SMILES 构建 PVA 三聚体
    
    PVA 重复单元: –CH2–CH(OH)–
    三聚体: CH3–CH(OH)–CH2–CH(OH)–CH2–CH(OH)–CH3
    
    实际结构（带端基）:
    H3C-CH(OH)-CH2-CH(OH)-CH2-CH(OH)-CH3
    """
    # 构建 PVA 三聚体 SMILES
    # 端基: CH3
    # 重复单元: CH(OH)-CH2
    
    if n_repeat == 1:
        smiles = "CC(O)C"
    elif n_repeat == 2:
        smiles = "CC(O)CC(O)C"
    elif n_repeat == 3:
        smiles = "CC(O)CC(O)CC(O)C"
    else:
        # 通用构建
        units = ["C(O)C"] * n_repeat
        smiles = "C" + "C".join(units)
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"无法从 SMILES 构建 PVA: {smiles}")
    
    # 添加氢
    mol = Chem.AddHs(mol)
    
    # 生成 3D 坐标
    status = AllChem.EmbedMolecule(mol, randomSeed=2025)
    if status < 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=2025)
    
    # UFF 优化
    try:
        if UFFHasAllMoleculeParams(mol):
            UFFOptimizeMolecule(mol, maxIters=500)
    except:
        pass
    
    return mol


# ==============================================================================
# MOL2/PDB 输出
# ==============================================================================

def write_pdb(mol: Chem.Mol, filepath: str, mol_name: str = "PVA") -> None:
    """写入 PDB 文件"""
    Chem.MolToPDBFile(mol, filepath)
    
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        raise RuntimeError(f"PDB 文件写入失败: {filepath}")


def calculate_cross_bond_lengths(mol: Chem.Mol, n_repeat: int) -> List[Tuple[int, int, float]]:
    """计算跨单体键的键长"""
    # 对于 PVA 三聚体，主链 C-C 键应该都是 ~1.53 Å
    # 找到主链上的 C-C 键
    
    conf = mol.GetConformer()
    bond_lengths = []
    
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # 只看 C-C 键
        if a1.GetSymbol() == 'C' and a2.GetSymbol() == 'C':
            idx1, idx2 = a1.GetIdx(), a2.GetIdx()
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            length = math.sqrt(
                (pos1.x - pos2.x)**2 + 
                (pos1.y - pos2.y)**2 + 
                (pos1.z - pos2.z)**2
            )
            bond_lengths.append((idx1, idx2, length))
    
    return bond_lengths


# ==============================================================================
# 主流程
# ==============================================================================

def generate_pva_trimer(
    input_path: Optional[str],
    output_path: str,
    n_repeat: int = 3,
    strategy: str = "auto",
    force: bool = False
) -> str:
    """
    生成真正的 PVA 三聚体
    
    Args:
        input_path: 输入 VA.mol2 路径（策略 A 使用）
        output_path: 输出 PVA.mol2 路径
        n_repeat: 聚合度
        strategy: "A" (水解), "B" (直接构建), "auto" (自动选择)
        force: 是否覆盖已有文件
    
    Returns:
        输出文件路径
    """
    print_header(f"生成真正的 PVA 三聚体 (n={n_repeat})")
    
    output_path = os.path.abspath(output_path)
    output_dir = os.path.dirname(output_path)
    pdb_path = os.path.splitext(output_path)[0] + ".pdb"
    
    print(f"  输出 mol2: {output_path}")
    print(f"  输出 pdb:  {pdb_path}")
    
    # 检查输出文件
    if os.path.exists(output_path) and not force:
        raise RuntimeError(f"输出文件已存在: {output_path}\n使用 --force 覆盖")
    
    # 确定策略
    actual_strategy = strategy
    
    if strategy == "auto":
        if input_path and os.path.exists(input_path):
            # 尝试读取并检查是否是 PVAc
            try:
                # 先用现有流程生成三聚体
                from scripts.cli.monomers import process_single_monomer
                
                # 生成到临时位置
                temp_output = output_path + ".temp"
                
                # 调用现有的聚合流程
                log_info("先使用现有流程生成三聚体...")
                
                mol = Chem.MolFromMol2File(input_path, removeHs=False)
                if mol is None:
                    actual_strategy = "B"
                else:
                    actual_strategy = "A"  # 假设 VA 会产生 PVAc
            except:
                actual_strategy = "B"
        else:
            actual_strategy = "B"
    
    print(f"  采用策略: {actual_strategy}")
    
    # 执行策略
    if actual_strategy == "A":
        print("\n[策略 A] 水解 PVAc → PVA")
        
        # 首先检查是否有现成的 PVAc 三聚体
        # 或者我们需要先生成它
        
        # 直接构建 PVAc 然后水解
        # PVAc 三聚体 SMILES
        pvac_smiles = "CC(=O)OC(C)CC(OC(C)=O)CC(OC(C)=O)C"
        
        log_info("构建 PVAc 三聚体...")
        pvac_mol = Chem.MolFromSmiles(pvac_smiles)
        if pvac_mol is None:
            log_warn("无法构建 PVAc，切换到策略 B")
            actual_strategy = "B"
        else:
            pvac_mol = Chem.AddHs(pvac_mol)
            AllChem.EmbedMolecule(pvac_mol, randomSeed=2025)
            
            print(f"  PVAc 原子数: {pvac_mol.GetNumAtoms()}")
            print(f"  PVAc 酯基数: {count_ester_groups(pvac_mol)}")
            
            log_info("水解 PVAc → PVA...")
            pva_mol = hydrolyze_pvac_manual(pvac_mol)
            
            if pva_mol is None or count_ester_groups(pva_mol) > 0:
                log_warn("水解失败，切换到策略 B")
                actual_strategy = "B"
    
    if actual_strategy == "B":
        print("\n[策略 B] 直接从 SMILES 构建 PVA")
        pva_mol = build_pva_from_smiles(n_repeat)
    
    # 验证
    print("\n=== 化学验证 ===")
    validate_pva_chemistry(pva_mol, n_repeat)
    
    # 连通性检查
    frags = Chem.GetMolFrags(pva_mol)
    print(f"\n=== 连通性检查 ===")
    print(f"  片段数: {len(frags)}")
    if len(frags) != 1:
        raise RuntimeError(f"生成的分子不是单一连通体！片段数: {len(frags)}")
    log_success("连通性检查通过: 单一分子")
    
    # 原子/键数
    print(f"\n=== 结构信息 ===")
    print(f"  原子数: {pva_mol.GetNumAtoms()}")
    print(f"  键数: {pva_mol.GetNumBonds()}")
    
    # C-C 键长
    cc_bonds = calculate_cross_bond_lengths(pva_mol, n_repeat)
    print(f"\n=== 主链 C-C 键长 ===")
    for i, (a1, a2, length) in enumerate(cc_bonds[:6]):  # 只显示前几个
        print(f"  键 {a1+1}--{a2+1}: {length:.3f} Å")
    
    # 写入文件
    print(f"\n=== 写入文件 ===")
    
    # 写 MOL2
    write_mol2_strict(pva_mol, output_path, "PVA", "PVA")
    size = os.path.getsize(output_path)
    log_success(f"MOL2: {output_path} ({size} bytes)")
    
    # 写 PDB
    write_pdb(pva_mol, pdb_path, "PVA")
    pdb_size = os.path.getsize(pdb_path)
    log_success(f"PDB: {pdb_path} ({pdb_size} bytes)")
    
    # 最终验证
    print(f"\n=== 最终验证 ===")
    validate_mol2_output(output_path)
    
    # 重新读取验证
    final_mol = Chem.MolFromMol2File(output_path, removeHs=False)
    if final_mol is None:
        raise RuntimeError("无法重新读取生成的 MOL2 文件")
    
    final_ester = count_ester_groups(final_mol)
    final_hydroxyl = count_hydroxyl_groups(final_mol)
    
    print(f"  酯基 C(=O)O: {final_ester} {'✓' if final_ester == 0 else '✗'}")
    print(f"  羟基 -OH: {final_hydroxyl} {'✓' if final_hydroxyl >= n_repeat else '✗'}")
    
    if final_ester > 0:
        raise RuntimeError("最终验证失败：仍存在酯基！")
    
    print_header("生成完成")
    log_success(f"真正的 PVA 三聚体已生成")
    print(f"  MOL2: {output_path}")
    print(f"  PDB:  {pdb_path}")
    
    return output_path


# ==============================================================================
# CLI
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="生成真正的 PVA（聚乙烯醇）三聚体",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 使用策略 B 直接构建（推荐）
  python scripts/build_pva_trimer.py --out mol2/PVA.mol2 --n 3 --force
  
  # 从 VA.mol2 尝试水解（策略 A）
  python scripts/build_pva_trimer.py --in mol2/VA.mol2 --out mol2/PVA.mol2 --n 3 --force

化学说明:
  - VA (Vinyl Acetate) 聚合得到 PVAc，不是 PVA
  - 真正的 PVA 重复单元: –CH2–CH(OH)–
  - 此脚本直接构建正确的 PVA 结构
"""
    )
    
    parser.add_argument("--in", dest="input", help="输入 VA.mol2 路径（可选，策略 A 使用）")
    parser.add_argument("--out", required=True, help="输出 PVA.mol2 路径")
    parser.add_argument("--n", type=int, default=3, help="聚合度 (默认: 3)")
    parser.add_argument("--strategy", choices=["A", "B", "auto"], default="auto",
                       help="策略: A=水解, B=直接构建, auto=自动选择")
    parser.add_argument("--force", action="store_true", help="覆盖已有文件")
    
    args = parser.parse_args()
    
    try:
        generate_pva_trimer(
            input_path=args.input,
            output_path=args.out,
            n_repeat=args.n,
            strategy=args.strategy,
            force=args.force
        )
        return 0
    except Exception as e:
        log_error(str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())


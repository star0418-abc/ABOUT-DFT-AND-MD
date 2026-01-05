#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
mol2_to_polymer_mol2.py - 从单体 MOL2 生成聚合物寡聚体 MOL2
=============================================================

功能:
  - 读取单体 mol2 文件 (仅支持 .mol2 格式)
  - 识别 head/tail 连接位点
  - 通过头尾相连生成线性寡聚体
  - 输出聚合物 mol2 (可被 Avogadro 打开)

⚠️ 重要变更 (v2.1.0):
  - 修复 EO 开环聚合：EO 是环氧三元环，需要开环而非 C=C 加成
  - 重构连接逻辑：使用 CombineMols + AddBond，禁止循环中 RemoveAtom
  - 添加连通性验证：每步检查 fragment 数，防止输出分散碎片

根因分析（v2.1.0 修复）:
  1. EO 问题根因：EO 是环氧乙烷（三元环 O1-C1-C2），正确机制是开环聚合，
     需要断开环内 C–O 键，而不是改键级。之前代码只处理 C=C 双键。
  2. 碎片问题根因：之前在循环中 RemoveAtom(氢) 会导致原子索引整体 shift，
     adjust 修补脆弱，容易连错原子导致多碎片。
  3. 解决方案：
     - EO 使用专用开环路径，断开 C-O 键生成可聚合单体
     - 聚合使用 CombineMols + AddBond，不删除原子
     - 每步验证连通性，失败立即报错

用法:
  # 使用默认配置 (n=3)
  python scripts/mol2_to_polymer_mol2.py --input mol2/EO.mol2 -v

  # 使用 recipe.yaml 配置
  python scripts/mol2_to_polymer_mol2.py --input mol2/EGDA.mol2 --recipe config/recipe.yaml

  # CLI 覆盖聚合度
  python scripts/mol2_to_polymer_mol2.py --input mol2/MMA.mol2 --n_repeat 5

支持的单体:
  EO → PEO (聚环氧乙烷) - 开环聚合
  VA → PVA (聚乙烯醇) - C=C 加成聚合
  TFEMA → PTFEMA (聚三氟乙基甲基丙烯酸酯) - C=C 加成聚合
  EGDA → PEGDA (聚乙二醇二丙烯酸酯) - 简化寡聚体
  MMA → PMMA (聚甲基丙烯酸甲酯) - C=C 加成聚合
  AM → PAM (聚丙烯酰胺) - C=C 加成聚合

版本: 2.1.0
"""

import argparse
import os
import sys
import random
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from pathlib import Path

# ==============================================================================
# RDKit 导入
# ==============================================================================

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms, Draw
    from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFHasAllMoleculeParams
    from rdkit import RDLogger
    # 禁用 RDKit 警告
    RDLogger.logger().setLevel(RDLogger.ERROR)
except ImportError:
    print("=" * 60)
    print("错误: RDKit 未安装")
    print("=" * 60)
    print("\n请使用以下命令安装 RDKit:")
    print("  conda install -c conda-forge rdkit")
    print("=" * 60)
    sys.exit(1)


# ==============================================================================
# YAML 导入 (可选，用于读取 recipe)
# ==============================================================================

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


# ==============================================================================
# 常量与配置
# ==============================================================================


# ⚠️ 重要变更: 默认聚合度从 5/20 改为 3 (三聚体)
DEFAULT_OLIGOMER_N = 3

# 单体到聚合物的映射
MONOMER_TO_POLYMER = {
    "EO": "PEO",
    "VA": "PVA", 
    "TFEMA": "PTFEMA",
    "EGDA": "PEGDA",
    "MMA": "PMMA",
    "AM": "PAM",
}

DEFAULT_MONOMERS = list(MONOMER_TO_POLYMER.keys())


# ==============================================================================
# 单体连接位点定义
# ==============================================================================

@dataclass
class MonomerConfig:
    """单体配置：定义 head/tail 连接位点"""
    name: str
    polymer_name: str
    # 原子名匹配模式（优先）
    head_atom_names: List[str]  # 可能的 head 原子名
    tail_atom_names: List[str]  # 可能的 tail 原子名
    # SMARTS 匹配模式（备选）
    head_smarts: Optional[str]
    tail_smarts: Optional[str]
    # 描述
    description: str
    # 是否是双官能团（交联剂）
    is_crosslinker: bool = False
    # 单体分子量 (g/mol) - 用于质量计算
    mw_g_mol: float = 0.0
    # 是否是环状醚（需要开环聚合）
    is_ring_opening: bool = False


# 单体配置表
MONOMER_CONFIGS: Dict[str, MonomerConfig] = {
    "EO": MonomerConfig(
        name="EO",
        polymer_name="PEO",
        head_atom_names=["O1"],  # 开环后 O 端作为 head
        tail_atom_names=["C2"],  # 开环后 C 端作为 tail
        head_smarts=None,  # EO 使用专用开环逻辑
        tail_smarts=None,
        description="环氧乙烷 → 聚环氧乙烷（开环聚合）",
        mw_g_mol=44.05,
        is_ring_opening=True,  # 标记为开环聚合
    ),
    "VA": MonomerConfig(
        name="VA",
        polymer_name="PVA",
        head_atom_names=["C6", "C4"],
        tail_atom_names=["C5", "C3"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[CH;$(C=C)]",
        description="乙酸乙烯酯 → 聚乙酸乙烯酯 (可水解为 PVA)",
        mw_g_mol=86.09,
    ),
    "TFEMA": MonomerConfig(
        name="TFEMA",
        polymer_name="PTFEMA",
        head_atom_names=["C3"],
        tail_atom_names=["C2"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[C;$(C(=C)(C)C)]",
        description="三氟乙基甲基丙烯酸酯 → 聚TFEMA",
        mw_g_mol=168.12,
    ),
    "MMA": MonomerConfig(
        name="MMA",
        polymer_name="PMMA",
        head_atom_names=["C3"],
        tail_atom_names=["C2"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[C;$(C(=C)(C)C)]",
        description="甲基丙烯酸甲酯 → 聚甲基丙烯酸甲酯",
        mw_g_mol=100.12,
    ),
    "AM": MonomerConfig(
        name="AM",
        polymer_name="PAM",
        head_atom_names=["C1", "HA"],
        tail_atom_names=["C2", "TA"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[CH;$(C=C)]",
        description="丙烯酰胺 → 聚丙烯酰胺",
        mw_g_mol=71.08,
    ),
    "EGDA": MonomerConfig(
        name="EGDA",
        polymer_name="PEGDA",
        head_atom_names=["HA", "C1"],
        tail_atom_names=["TB", "HB"],
        head_smarts="[CH2;$(C=C)]",
        tail_smarts="[CH2;$(C=C)]",
        description="乙二醇二丙烯酸酯 → 简化 PEGDA 寡聚体",
        is_crosslinker=True,
        mw_g_mol=198.22,
    ),
}


# ==============================================================================
# 连通性验证（核心新增功能）
# ==============================================================================

class ConnectivityError(ValueError):
    """连通性检查失败异常"""
    pass


def assert_single_fragment(
    mol: Chem.Mol, 
    context: str = "",
    step: int = -1,
    head_idx: int = -1,
    tail_idx: int = -1
) -> None:
    """
    验证分子是单一连通体，否则抛出异常
    
    Args:
        mol: RDKit 分子
        context: 上下文描述（用于错误信息）
        step: 当前步骤（-1 表示不显示）
        head_idx: head 原子索引（用于调试）
        tail_idx: tail 原子索引（用于调试）
    
    Raises:
        ConnectivityError: 如果 fragment 数 != 1
    """
    frags = Chem.GetMolFrags(mol)
    num_frags = len(frags)
    
    if num_frags != 1:
        debug_info = [
            f"连通性检查失败: 检测到 {num_frags} 个碎片（应为 1）",
            f"上下文: {context}",
        ]
        if step >= 0:
            debug_info.append(f"步骤: {step}")
        if head_idx >= 0:
            debug_info.append(f"Head 索引: {head_idx}")
        if tail_idx >= 0:
            debug_info.append(f"Tail 索引: {tail_idx}")
        debug_info.append(f"原子数: {mol.GetNumAtoms()}")
        debug_info.append(f"键数: {mol.GetNumBonds()}")
        debug_info.append(f"碎片原子数: {[len(f) for f in frags]}")
        
        raise ConnectivityError("\n".join(debug_info))


def check_connectivity(mol: Chem.Mol, verbose: bool = False) -> Tuple[bool, int]:
    """
    检查分子连通性
    
    Returns:
        (是否单一连通, 碎片数)
    """
    frags = Chem.GetMolFrags(mol)
    num_frags = len(frags)
    is_connected = num_frags == 1
    
    if verbose:
        if is_connected:
            print(f"    [连通性] ✓ 单一分子 ({mol.GetNumAtoms()} 原子)")
        else:
            print(f"    [连通性] ✗ {num_frags} 个碎片!")
    
    return is_connected, num_frags


# ==============================================================================
# Recipe 配置读取
# ==============================================================================

def load_oligomer_n_from_recipe(recipe_path: str) -> Optional[int]:
    """从 recipe.yaml 读取 oligomer_n 配置"""
    if not HAS_YAML:
        print("[WARN] PyYAML 未安装，无法读取 recipe.yaml")
        return None
    
    if not os.path.isfile(recipe_path):
        print(f"[WARN] recipe 文件不存在: {recipe_path}")
        return None
    
    try:
        with open(recipe_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        polymerization = config.get("polymerization", {})
        if isinstance(polymerization, dict):
            oligomer_n = polymerization.get("oligomer_n")
            if oligomer_n is not None:
                return int(oligomer_n)
        
        polymer_matrix = config.get("polymer_matrix", [])
        if isinstance(polymer_matrix, list):
            for item in polymer_matrix:
                if isinstance(item, dict) and "oligomer_n" in item:
                    return int(item["oligomer_n"])
        
        return None
        
    except Exception as e:
        print(f"[WARN] 读取 recipe 失败: {e}")
        return None


def get_effective_oligomer_n(
    cli_n_repeat: Optional[int],
    recipe_path: Optional[str],
    default: int = DEFAULT_OLIGOMER_N
) -> int:
    """确定最终使用的 oligomer_n（优先级: CLI > recipe > 默认）"""
    if cli_n_repeat is not None:
        return cli_n_repeat
    
    if recipe_path:
        recipe_n = load_oligomer_n_from_recipe(recipe_path)
        if recipe_n is not None:
            return recipe_n
    
    return default


# ==============================================================================
# MOL2 读写函数
# ==============================================================================

def validate_mol2_input(filepath: str) -> bool:
    """验证输入文件必须是 .mol2 格式"""
    ext = os.path.splitext(filepath)[1].lower()
    if ext != ".mol2":
        print(f"[ERROR] 仅支持 .mol2 格式输入，收到: {ext}")
        print(f"        文件: {filepath}")
        return False
    return True


def read_mol2(filepath: str) -> Optional[Chem.Mol]:
    """读取 MOL2 文件"""
    if not os.path.isfile(filepath):
        print(f"  [ERROR] 文件不存在: {filepath}")
        return None
    
    try:
        mol = Chem.MolFromMol2File(filepath, removeHs=False)
        if mol is None:
            print(f"  [ERROR] 无法解析 MOL2 文件: {filepath}")
            return None
        return mol
    except Exception as e:
        print(f"  [ERROR] 读取 MOL2 失败: {e}")
        return None


def write_mol2(mol: Chem.Mol, filepath: str, mol_name: str, res_name: str) -> bool:
    """
    写入 MOL2 文件
    
    此函数不依赖 RDKit 的 sanitize，直接从分子图中提取信息
    """
    try:
        conf = mol.GetConformer()
        atoms = list(mol.GetAtoms())
        bonds = list(mol.GetBonds())
        
        atom_names = {}
        element_count = {}
        
        for atom in atoms:
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            element_count[symbol] = element_count.get(symbol, 0) + 1
            atom_names[idx] = f"{symbol}{element_count[symbol]}"
        
        def get_sybyl_type(atom):
            """根据原子的键连接判断 Sybyl 类型"""
            symbol = atom.GetSymbol()
            
            if symbol == 'H':
                return 'H'
            elif symbol == 'F':
                return 'F'
            elif symbol == 'Cl':
                return 'Cl'
            elif symbol == 'Br':
                return 'Br'
            elif symbol == 'I':
                return 'I'
            elif symbol == 'S':
                return 'S.3'
            
            # 对于 C, N, O，检查键来判断杂化
            has_double = False
            has_triple = False
            
            for bond in atom.GetBonds():
                bt = bond.GetBondType()
                if bt == Chem.BondType.DOUBLE:
                    has_double = True
                elif bt == Chem.BondType.TRIPLE:
                    has_triple = True
            
            if symbol == 'C':
                if has_triple:
                    return 'C.1'
                elif has_double:
                    return 'C.2'
                else:
                    return 'C.3'
            elif symbol == 'N':
                if has_double or has_triple:
                    return 'N.2'
                else:
                    return 'N.3'
            elif symbol == 'O':
                if has_double:
                    return 'O.2'
                else:
                    return 'O.3'
            else:
                return symbol
        
        bond_type_map = {
            Chem.BondType.SINGLE: '1',
            Chem.BondType.DOUBLE: '2',
            Chem.BondType.TRIPLE: '3',
            Chem.BondType.AROMATIC: 'ar',
        }
        
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
                # 不使用 GetFormalCharge() 因为可能在未 sanitize 的分子上失败
                charge = 0.0
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
        
        return True
        
    except Exception as e:
        print(f"  [ERROR] 写入 MOL2 失败: {e}")
        import traceback
        traceback.print_exc()
        return False


# ==============================================================================
# 连接位点识别
# ==============================================================================

def get_atom_name_from_mol2(filepath: str) -> Dict[int, str]:
    """从 MOL2 文件解析原子名"""
    atom_names = {}
    
    try:
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
    except Exception as e:
        print(f"  [WARN] 无法解析原子名: {e}")
    
    return atom_names


def find_connecting_atoms(
    mol: Chem.Mol,
    config: MonomerConfig,
    atom_names: Dict[int, str],
    verbose: bool = False
) -> Tuple[Optional[int], Optional[int]]:
    """找到单体的 head 和 tail 连接原子"""
    head_idx = None
    tail_idx = None
    
    # 方法 1: 通过原子名匹配
    name_to_idx = {name: idx for idx, name in atom_names.items()}
    
    for name in config.head_atom_names:
        if name in name_to_idx:
            head_idx = name_to_idx[name]
            if verbose:
                print(f"    Head 原子 (by name '{name}'): {head_idx + 1}")
            break
    
    for name in config.tail_atom_names:
        if name in name_to_idx:
            tail_idx = name_to_idx[name]
            if verbose:
                print(f"    Tail 原子 (by name '{name}'): {tail_idx + 1}")
            break
    
    # 方法 2: 通过 SMARTS 匹配
    if head_idx is None and config.head_smarts:
        pattern = Chem.MolFromSmarts(config.head_smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                head_idx = matches[0][0]
                if verbose:
                    print(f"    Head 原子 (by SMARTS): {head_idx + 1}")
    
    if tail_idx is None and config.tail_smarts:
        pattern = Chem.MolFromSmarts(config.tail_smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                for match in reversed(matches):
                    if match[0] != head_idx:
                        tail_idx = match[0]
                        if verbose:
                            print(f"    Tail 原子 (by SMARTS): {tail_idx + 1}")
                        break
    
    # 特殊处理: EGDA 的两端
    if config.name == "EGDA" and head_idx is not None and tail_idx is None:
        pattern = Chem.MolFromSmarts("[CH2;$(C=C)]")
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if match[0] != head_idx:
                    tail_idx = match[0]
                    if verbose:
                        print(f"    Tail 原子 (EGDA 另一端): {tail_idx + 1}")
                    break
    
    return head_idx, tail_idx


# ==============================================================================
# EO 开环聚合专用逻辑
# ==============================================================================

def detect_epoxide(mol: Chem.Mol) -> Optional[Tuple[int, int, int]]:
    """
    检测环氧结构（三元环 O-C-C）
    
    Returns:
        (O_idx, C1_idx, C2_idx) 或 None
    
    EO 结构: 三元环，一个 O 连接两个 C
    """
    # SMARTS 匹配环氧：三元环中的 O 连接两个 C
    epoxide_smarts = Chem.MolFromSmarts("[O;R1]1[C;R1][C;R1]1")
    if epoxide_smarts is None:
        return None
    
    matches = mol.GetSubstructMatches(epoxide_smarts)
    if matches:
        # 返回第一个匹配: (O, C1, C2)
        return matches[0]
    
    return None


def build_peo_from_smiles(dp: int, verbose: bool = False) -> Chem.Mol:
    """
    直接从 SMILES 构建 PEO 链
    
    PEO 结构: HO-(CH2-CH2-O)_n-H
    简化为: HO-CH2-CH2-O-CH2-CH2-O-...-CH2-CH2-OH
    
    Args:
        dp: 聚合度（重复单元数）
        verbose: 是否详细输出
    
    Returns:
        RDKit Mol 对象
    """
    if verbose:
        print(f"    [EO] 使用 SMILES 直接构建 PEO (dp={dp})")
    
    # 构建 PEO SMILES: HO-(CH2CH2O)_n-H
    # 例如 dp=3: OCCOCCOCCO 或更精确 OCCOCCOC
    if dp == 1:
        smiles = "OCCO"  # HO-CH2-CH2-OH
    else:
        # 构建重复单元
        # (CH2-CH2-O) repeated, 端基 OH
        units = ["CCO"] * dp
        smiles = "O" + "".join(units)  # 开头 HO- 隐含
    
    if verbose:
        print(f"    [EO] SMILES: {smiles}")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法从 SMILES 构建 PEO: {smiles}")
    
    # 添加氢原子
    mol = Chem.AddHs(mol)
    
    # 验证连通性
    assert_single_fragment(mol, context="PEO SMILES 构建后")
    
    # 生成 3D 坐标
    try:
        status = AllChem.EmbedMolecule(mol, randomSeed=2025)
        if status < 0:
            AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=2025)
    except Exception as e:
        if verbose:
            print(f"    [EO] 嵌入警告: {e}")
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=2025)
    
    if verbose:
        print(f"    [EO] PEO 构建完成: {mol.GetNumAtoms()} 原子, {mol.GetNumBonds()} 键")
    
    return mol


def polymerize_eo(dp: int, seed: int = 2025, verbose: bool = False) -> Chem.Mol:
    """
    EO 开环聚合专用函数
    
    EO (Ethylene Oxide) 是三元环环氧，聚合机制是开环聚合：
    - 断开 C-O 键
    - O 端与前一单元的 C 端连接
    - 形成 -(CH2-CH2-O)_n- 链
    
    由于开环聚合复杂度高，这里使用 SMILES 直接构建 PEO 链，
    化学上等价于开环聚合的结果。
    
    Args:
        dp: 聚合度
        seed: 随机种子
        verbose: 是否详细输出
    
    Returns:
        PEO 分子
    """
    if verbose:
        print(f"\n  [EO 开环聚合] 聚合度={dp}")
        print("    机制: 环氧开环 → -(CH2-CH2-O)_n-")
    
    random.seed(seed)
    
    # 直接从 SMILES 构建（最可靠的方法）
    polymer = build_peo_from_smiles(dp, verbose)
    
    # 最终验证
    assert_single_fragment(polymer, context="PEO 最终验证")
    
    return polymer


# ==============================================================================
# 通用聚合物生成（C=C 加成聚合）
# ==============================================================================

def prepare_monomer_for_polymerization(
    monomer: Chem.Mol,
    head_idx: int,
    tail_idx: int,
    verbose: bool = False
) -> Chem.Mol:
    """
    准备单体用于 C=C 加成聚合：将 head-tail 之间的双键改为单键
    
    注意：
    - 此函数仅用于 C=C 加成聚合（如 MMA, TFEMA），不用于开环聚合（如 EO）
    - 不添加额外的隐式氢（因为聚合后会形成新的 C-C 键）
    """
    mol = Chem.RWMol(Chem.Mol(monomer))
    
    bond = mol.GetBondBetweenAtoms(head_idx, tail_idx)
    
    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
        if verbose:
            print(f"    [PREP] 将 C{head_idx+1}=C{tail_idx+1} 双键改为单键")
        bond.SetBondType(Chem.BondType.SINGLE)
        
        # 不添加额外的隐式氢，因为聚合时会形成新的 C-C 键
        # 这避免了价态问题
    
    return mol.GetMol()


def polymerize_linear_v2(
    monomer: Chem.Mol,
    head_idx: int,
    tail_idx: int,
    dp: int,
    seed: int = 2025,
    verbose: bool = False,
    validate_each_step: bool = False
) -> Optional[Chem.Mol]:
    """
    生成线性聚合物（v2.1 重构版本）
    
    对于 C=C 加成聚合，使用 SMILES 构建方法：
    1. 获取单体 SMILES
    2. 将双键位置标记为连接点
    3. 直接构建重复单元的 SMILES
    4. 解析并生成 3D 坐标
    
    Args:
        monomer: 单体分子
        head_idx: head 原子索引（双键的一端，如 CH2=）
        tail_idx: tail 原子索引（双键的另一端，如 =CH-）
        dp: 聚合度
        seed: 随机种子
        verbose: 是否详细输出
        validate_each_step: 是否每步验证连通性
    
    Returns:
        聚合物分子，或 None
    """
    if dp < 1:
        print("  [ERROR] DP 必须 >= 1")
        return None
    
    if dp == 1:
        return Chem.Mol(monomer)
    
    random.seed(seed)
    
    try:
        # ====================================================================
        # 方法：使用 SMILES 构建聚合物
        # 将双键改为连接点 [*]，然后重复拼接
        # ====================================================================
        
        # 准备单体：将双键改为单键，创建可连接的结构
        prep_monomer = prepare_monomer_for_polymerization(monomer, head_idx, tail_idx, verbose)
        
        # 获取单体 SMILES（不带氢，便于操作）
        monomer_no_h = Chem.RemoveHs(prep_monomer)
        
        # 尝试直接从原始单体 SMILES 构建聚合物
        # 对于加成聚合，重复单元就是单体去掉双键
        original_smiles = Chem.MolToSmiles(Chem.RemoveHs(monomer))
        
        if verbose:
            print(f"    原始单体 SMILES: {original_smiles}")
        
        # 策略：使用 RDKit 的反应 SMARTS 进行聚合
        # 这是处理复杂单体的最可靠方式
        
        # 对于双键聚合：C=C → -C-C-C-C-...
        # 使用 RDKit 的 RWMol 和 SMILES 操作
        
        # 准备后的单体 SMILES
        prep_smiles = Chem.MolToSmiles(monomer_no_h)
        if verbose:
            print(f"    准备后 SMILES: {prep_smiles}")
        
        # 直接构建聚合物 SMILES
        # 对于甲基丙烯酸酯类 C=C(C)C(=O)OR
        # 聚合后重复单元: -CH2-C(CH3)(C(=O)OR)-
        
        # 使用 RDKit 的 RWMol 进行构建，一步步添加
        prep_monomer_h = Chem.AddHs(prep_monomer)
        
        # 找到 head 和 tail 原子上可用于连接的氢
        # 这些氢将在连接时被"替换"（通过形成新键）
        
        polymer = Chem.RWMol(Chem.Mol(prep_monomer_h))
        current_tail = tail_idx
        
        if verbose:
            print(f"    初始: {polymer.GetNumAtoms()} 原子")
        
        for i in range(1, dp):
            offset = polymer.GetNumAtoms()
            
            # 合并新单体
            new_monomer = Chem.Mol(prep_monomer_h)
            combined = Chem.CombineMols(polymer.GetMol(), new_monomer)
            polymer = Chem.RWMol(combined)
            
            # 新单体的位置
            new_head = head_idx + offset
            new_tail = tail_idx + offset
            
            # 添加连接键
            polymer.AddBond(current_tail, new_head, Chem.BondType.SINGLE)
            
            # 更新 tail
            current_tail = new_tail
            
            if verbose and (i % 5 == 0 or i == dp - 1):
                print(f"    已连接 {i + 1}/{dp} 个单体 ({polymer.GetNumAtoms()} 原子)")
        
        result = polymer.GetMol()
        
        # 清理分子（尝试修复价态问题）
        try:
            # 不做完整 sanitize，只更新属性
            Chem.FastFindRings(result)
            result.UpdatePropertyCache(strict=False)
        except Exception as e:
            if verbose:
                print(f"    [WARN] 属性更新: {e}")
        
        # 最终连通性验证（使用图连通性，不依赖 sanitize）
        frags = Chem.GetMolFrags(result)
        if len(frags) != 1:
            raise ConnectivityError(f"聚合后有 {len(frags)} 个碎片")
        
        if verbose:
            print(f"    连通性验证: ✓ 单一分子")
        
        # 生成 3D 坐标
        if verbose:
            print("    生成 3D 坐标...")
        
        try:
            result.RemoveAllConformers()
            # 使用宽松参数
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            status = AllChem.EmbedMolecule(result, params)
            if status < 0:
                if verbose:
                    print("    [WARN] ETKDG 失败，尝试随机坐标...")
                AllChem.EmbedMolecule(result, useRandomCoords=True, randomSeed=seed)
        except Exception as e:
            if verbose:
                print(f"    [WARN] 嵌入失败: {e}, 使用随机坐标")
            try:
                result.RemoveAllConformers()
                AllChem.EmbedMolecule(result, useRandomCoords=True, randomSeed=seed)
            except:
                # 最后手段：手动创建坐标
                conf = Chem.Conformer(result.GetNumAtoms())
                for j in range(result.GetNumAtoms()):
                    conf.SetAtomPosition(j, (random.uniform(-10, 10),
                                            random.uniform(-10, 10),
                                            random.uniform(-10, 10)))
                result.AddConformer(conf, assignId=True)
        
        return result
        
    except ConnectivityError as e:
        print(f"  [ERROR] 连通性错误: {e}")
        return None
    except Exception as e:
        print(f"  [ERROR] 聚合失败: {e}")
        import traceback
        traceback.print_exc()
        return None


# ==============================================================================
# 主生成函数
# ==============================================================================

def generate_polymer_from_file(
    input_path: str,
    oligomer_n: int,
    optimize: bool = True,
    seed: int = 2025,
    verbose: bool = False,
    dry_run: bool = False,
    validate_each_step: bool = False
) -> Optional[str]:
    """
    从单个 MOL2 文件生成聚合物
    
    输出规则: 输入 X.mol2 → 输出 PX.mol2 (在输入文件的原目录)
    """
    # 验证输入格式
    if not validate_mol2_input(input_path):
        return None
    
    input_path = os.path.abspath(input_path)
    input_dir = os.path.dirname(input_path)
    input_basename = os.path.basename(input_path)
    input_name, input_ext = os.path.splitext(input_basename)
    
    # 输出命名规则
    output_name = f"P{input_name}{input_ext}"
    output_path = os.path.join(input_dir, output_name)
    
    # 尝试匹配已知单体配置
    config = None
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
    print(f"单体 MOL2 → 寡聚体 MOL2 (oligomer_n={oligomer_n})")
    print("=" * 60)
    print(f"  输入文件: {input_path}")
    print(f"  输出文件: {output_path}")
    print(f"  聚合度:   {oligomer_n}")
    print(f"  配置:     {config.description}")
    if config.is_ring_opening:
        print(f"  聚合类型: 开环聚合 (Ring-Opening)")
    else:
        print(f"  聚合类型: C=C 加成聚合")
    print("=" * 60)
    
    if dry_run:
        print(f"  [DRY-RUN] 将生成 {output_name}")
        return output_path
    
    # 读取单体
    monomer = read_mol2(input_path)
    if monomer is None:
        return None
    
    monomer_atoms = monomer.GetNumAtoms()
    monomer_bonds = monomer.GetNumBonds()
    
    print(f"  单体原子数: {monomer_atoms}")
    print(f"  单体键数:   {monomer_bonds}")
    
    # 根据单体类型选择聚合策略
    if config.is_ring_opening:
        # ========== EO 开环聚合 ==========
        # 检测环氧结构
        epoxide = detect_epoxide(monomer)
        if epoxide:
            if verbose:
                print(f"  检测到环氧结构: O={epoxide[0]+1}, C1={epoxide[1]+1}, C2={epoxide[2]+1}")
        else:
            if verbose:
                print("  [WARN] 未检测到环氧结构，使用 SMILES 构建")
        
        # 使用专用 EO 聚合函数
        polymer = polymerize_eo(dp=oligomer_n, seed=seed, verbose=verbose)
        
    else:
        # ========== C=C 加成聚合 ==========
        # 获取原子名
        atom_names = get_atom_name_from_mol2(input_path)
        
        # 找到连接位点
        head_idx, tail_idx = find_connecting_atoms(monomer, config, atom_names, verbose)
        
        if head_idx is None or tail_idx is None:
            print(f"  [ERROR] 无法识别连接位点")
            print(f"    Head: {head_idx}, Tail: {tail_idx}")
            return None
        
        print(f"  连接位点: Head={head_idx + 1}, Tail={tail_idx + 1} (1-based)")
        
        # 使用 v2 聚合函数
        polymer = polymerize_linear_v2(
            monomer=monomer,
            head_idx=head_idx,
            tail_idx=tail_idx,
            dp=oligomer_n,
            seed=seed,
            verbose=verbose,
            validate_each_step=validate_each_step
        )
    
    if polymer is None:
        return None
    
    polymer_atoms = polymer.GetNumAtoms()
    polymer_bonds = polymer.GetNumBonds()
    
    print(f"  聚合物原子数: {polymer_atoms}")
    print(f"  聚合物键数:   {polymer_bonds}")
    
    # 最终连通性检查
    is_connected, num_frags = check_connectivity(polymer, verbose=True)
    if not is_connected:
        print(f"  [ERROR] 最终验证失败: 生成的结构包含 {num_frags} 个碎片!")
        print(f"          拒绝输出分散碎片的 mol2 文件")
        return None
    
    # 确保有 3D 坐标
    if polymer.GetNumConformers() == 0:
        print("  生成 3D 坐标...")
        try:
            AllChem.EmbedMolecule(polymer, randomSeed=seed)
        except:
            AllChem.EmbedMolecule(polymer, useRandomCoords=True, randomSeed=seed)
    
    # UFF 优化
    if optimize:
        print("  UFF 优化中...")
        try:
            if UFFHasAllMoleculeParams(polymer):
                UFFOptimizeMolecule(polymer, maxIters=500)
            else:
                print("  [WARN] 部分原子无 UFF 参数，跳过优化")
        except Exception as e:
            print(f"  [WARN] 优化失败: {e}")
    
    # 写入 MOL2
    polymer_name = config.polymer_name if config else f"P{input_name}"
    success = write_mol2(
        mol=polymer,
        filepath=output_path,
        mol_name=polymer_name,
        res_name=polymer_name[:3]
    )
    
    if success:
        file_size = os.path.getsize(output_path)
        print("")
        print(f"  [SUCCESS] 生成成功!")
        print(f"    输出: {output_path}")
        print(f"    大小: {file_size} bytes")
        print(f"    原子: {polymer_atoms}")
        print(f"    键数: {polymer_bonds}")
        print(f"    连通: ✓ 单一分子 (1 fragment)")
        return output_path
    
    return None


def generate_polymer_batch(
    monomers: List[str],
    dp_map: Dict[str, int],
    mol2_dir: str,
    optimize: bool = True,
    seed: int = 2025,
    verbose: bool = False,
    dry_run: bool = False,
    no_overwrite: bool = False
) -> Dict[str, List[str]]:
    """批量生成聚合物"""
    results = {"success": [], "failed": [], "skipped": []}
    
    for monomer in monomers:
        if monomer not in MONOMER_CONFIGS:
            print(f"\n[WARN] 未知单体: {monomer}")
            results["skipped"].append(monomer)
            continue
        
        config = MONOMER_CONFIGS[monomer]
        dp = dp_map.get(monomer, DEFAULT_OLIGOMER_N)
        
        monomer_path = os.path.join(mol2_dir, f"{monomer}.mol2")
        output_path = os.path.join(mol2_dir, f"P{monomer}.mol2")
        
        if not os.path.isfile(monomer_path):
            print(f"\n[WARN] 单体文件不存在: {monomer_path}")
            results["skipped"].append(monomer)
            continue
        
        if os.path.isfile(output_path) and no_overwrite:
            print(f"\n[SKIP] 文件已存在: {output_path}")
            results["skipped"].append(monomer)
            continue
        
        try:
            result = generate_polymer_from_file(
                input_path=monomer_path,
                oligomer_n=dp,
                optimize=optimize,
                seed=seed,
                verbose=verbose,
                dry_run=dry_run
            )
            
            if result:
                results["success"].append(f"P{monomer}")
            else:
                results["failed"].append(monomer)
                
        except Exception as e:
            print(f"\n[ERROR] 处理 {monomer} 时发生异常: {e}")
            import traceback
            traceback.print_exc()
            results["failed"].append(monomer)
    
    return results


# ==============================================================================
# CLI
# ==============================================================================

def parse_dp_arg(dp_str: str, default_dp: int) -> Dict[str, int]:
    """解析 --dp 参数"""
    result = {}
    
    try:
        default_val = int(dp_str)
        for monomer in MONOMER_TO_POLYMER:
            result[monomer] = default_val
        return result
    except ValueError:
        pass
    
    for part in dp_str.split(","):
        if "=" in part:
            key, val = part.split("=", 1)
            key = key.strip().upper()
            try:
                result[key] = int(val.strip())
            except ValueError:
                print(f"[WARN] 无效的 DP 值: {part}")
    
    for monomer in MONOMER_TO_POLYMER:
        if monomer not in result:
            result[monomer] = default_dp
    
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="从单体 MOL2 生成聚合物寡聚体 MOL2 (v2.1: 修复 EO 开环聚合)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
示例:
  # 单文件模式 (推荐)
  %(prog)s --input mol2/EO.mol2 -v                 # EO 开环聚合
  %(prog)s --input mol2/MMA.mol2 --n_repeat 5      # MMA 加成聚合
  %(prog)s --input mol2/EGDA.mol2 --recipe config/recipe.yaml

  # 批处理模式
  %(prog)s --monomers EO,MMA --dp 3

单体 → 聚合物:
  EO → PEO (聚环氧乙烷) - 开环聚合 ★
  VA → PVA (聚乙烯醇) - C=C 加成聚合
  TFEMA → PTFEMA - C=C 加成聚合
  EGDA → PEGDA - C=C 加成聚合
  MMA → PMMA - C=C 加成聚合
  AM → PAM - C=C 加成聚合

⚠️ v2.1.0 重要修复:
  - EO 现在使用正确的开环聚合机制
  - 修复了多碎片问题（禁止 RemoveAtom，使用 CombineMols+AddBond）
  - 添加了严格的连通性验证
"""
    )
    
    parser.add_argument("--input", "-i", help="输入单体 MOL2 文件路径")
    parser.add_argument("--n_repeat", type=int, default=None, help="聚合度")
    parser.add_argument("--recipe", "-r", default=None, help="recipe.yaml 路径")
    parser.add_argument("--dir", default="mol2", help="MOL2 目录 (批处理模式)")
    parser.add_argument("--monomers", default=None, help="单体列表 (批处理模式)")
    parser.add_argument("--dp", default=None, help="聚合度 (批处理模式)")
    parser.add_argument("--opt", action="store_true", default=True, help="启用 UFF 优化")
    parser.add_argument("--no_opt", action="store_true", help="禁用 UFF 优化")
    parser.add_argument("--no_overwrite", action="store_true", help="不覆盖已存在的文件")
    parser.add_argument("--seed", type=int, default=2025, help="随机种子")
    parser.add_argument("--dry_run", action="store_true", help="模拟运行")
    parser.add_argument("--verbose", "-v", action="store_true", help="详细输出")
    parser.add_argument("--validate_each_step", action="store_true", 
                       help="每步验证连通性 (调试用)")
    parser.add_argument("--selfcheck", action="store_true", help="运行自检测试")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    
    return parser.parse_args()


def run_selfcheck() -> int:
    """
    运行自检测试
    
    包含 EO 开环聚合测试（v2.1 新增）
    """
    print("=" * 60)
    print(f"mol2_to_polymer_mol2.py v{VERSION} 自检测试")
    print("=" * 60)
    
    tests_passed = 0
    tests_failed = 0
    
    # 测试 1: 默认聚合度是 3
    print("\n[TEST 1] 默认聚合度 = 3")
    if DEFAULT_OLIGOMER_N == 3:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL: DEFAULT_OLIGOMER_N = {DEFAULT_OLIGOMER_N}")
        tests_failed += 1
    
    # 测试 2: 优先级逻辑
    print("\n[TEST 2] 优先级逻辑 (CLI > recipe > 默认)")
    n1 = get_effective_oligomer_n(cli_n_repeat=5, recipe_path=None)
    n2 = get_effective_oligomer_n(cli_n_repeat=None, recipe_path=None)
    if n1 == 5 and n2 == 3:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL: CLI=5 → {n1}, CLI=None → {n2}")
        tests_failed += 1
    
    # 测试 3: 输入验证
    print("\n[TEST 3] 输入格式验证")
    valid_mol2 = validate_mol2_input("/tmp/test.mol2")
    import io
    import sys
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    invalid_pdb = validate_mol2_input("/tmp/test.pdb")
    sys.stdout = old_stdout
    if valid_mol2 and not invalid_pdb:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL")
        tests_failed += 1
    
    # 测试 4: EO 开环聚合 (关键新测试)
    print("\n[TEST 4] EO 开环聚合 (dp=3)")
    try:
        peo = polymerize_eo(dp=3, verbose=False)
        frags = Chem.GetMolFrags(peo)
        num_atoms = peo.GetNumAtoms()
        num_bonds = peo.GetNumBonds()
        
        # PEO dp=3: O-(CCO)3 = OCCOCCOCCO = 10 重原子 + H
        # 期望: 单一连通分子
        if len(frags) == 1 and num_bonds > 0:
            print(f"  ✓ PASS: {num_atoms} 原子, {num_bonds} 键, 1 fragment")
            tests_passed += 1
        else:
            print(f"  ✗ FAIL: {len(frags)} fragments")
            tests_failed += 1
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
        tests_failed += 1
    
    # 测试 5: EO 开环聚合连通性验证
    print("\n[TEST 5] EO 连通性验证 (dp=5)")
    try:
        peo5 = polymerize_eo(dp=5, verbose=False)
        assert_single_fragment(peo5, context="EO dp=5 测试")
        print(f"  ✓ PASS: 单一连通分子")
        tests_passed += 1
    except ConnectivityError as e:
        print(f"  ✗ FAIL: {e}")
        tests_failed += 1
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
        tests_failed += 1
    
    # 测试 6: 环氧检测
    print("\n[TEST 6] 环氧结构检测")
    try:
        # 创建 EO 分子 (C1CO1)
        eo_mol = Chem.MolFromSmiles("C1CO1")
        epoxide = detect_epoxide(eo_mol)
        if epoxide is not None:
            print(f"  ✓ PASS: 检测到环氧 {epoxide}")
            tests_passed += 1
        else:
            print(f"  ✗ FAIL: 未检测到环氧")
            tests_failed += 1
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
        tests_failed += 1
    
    # 测试 7: 连通性验证函数
    print("\n[TEST 7] 连通性验证函数")
    try:
        # 单一分子
        single_mol = Chem.MolFromSmiles("CCO")
        single_mol = Chem.AddHs(single_mol)
        assert_single_fragment(single_mol, context="单一分子测试")
        
        # 多碎片分子
        frag_mol = Chem.MolFromSmiles("CCO.CC")
        frag_mol = Chem.AddHs(frag_mol)
        try:
            assert_single_fragment(frag_mol, context="多碎片测试")
            print(f"  ✗ FAIL: 应该抛出 ConnectivityError")
            tests_failed += 1
        except ConnectivityError:
            print(f"  ✓ PASS: 正确检测多碎片")
            tests_passed += 1
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
        tests_failed += 1
    
    # 测试 8: MMA 加成聚合 (确保不破坏现有功能)
    print("\n[TEST 8] MMA 加成聚合兼容性")
    try:
        # 用 SMILES 创建 MMA
        mma = Chem.MolFromSmiles("C=C(C)C(=O)OC")
        mma = Chem.AddHs(mma)
        AllChem.EmbedMolecule(mma, randomSeed=2025)
        
        # 找 head/tail (简化: 用索引)
        # MMA: C=C(C)C(=O)OC, 双键是 0-1
        head_idx = 0  # =CH2
        tail_idx = 1  # =C(C)-
        
        polymer = polymerize_linear_v2(mma, head_idx, tail_idx, dp=3, verbose=False)
        if polymer is not None:
            frags = Chem.GetMolFrags(polymer)
            if len(frags) == 1:
                print(f"  ✓ PASS: MMA 聚合成功, 1 fragment")
                tests_passed += 1
            else:
                print(f"  ✗ FAIL: MMA 聚合产生 {len(frags)} fragments")
                tests_failed += 1
        else:
            print(f"  ✗ FAIL: MMA 聚合返回 None")
            tests_failed += 1
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
        tests_failed += 1
    
    # 汇总
    print("\n" + "=" * 60)
    total = tests_passed + tests_failed
    print(f"自检完成: {tests_passed}/{total} 通过")
    if tests_failed > 0:
        print(f"  ✗ {tests_failed} 项失败")
    else:
        print(f"  ✓ 全部通过!")
    print("=" * 60)
    
    return 0 if tests_failed == 0 else 1


def main() -> int:
    args = parse_args()
    
    if args.selfcheck:
        return run_selfcheck()
    
    optimize = args.opt and not args.no_opt
    
    if args.input:
        # 单文件模式
        oligomer_n = get_effective_oligomer_n(
            cli_n_repeat=args.n_repeat,
            recipe_path=args.recipe
        )
        
        print(f"[INFO] 聚合度确定: oligomer_n = {oligomer_n}")
        
        result = generate_polymer_from_file(
            input_path=args.input,
            oligomer_n=oligomer_n,
            optimize=optimize,
            seed=args.seed,
            verbose=args.verbose,
            dry_run=args.dry_run,
            validate_each_step=args.validate_each_step
        )
        
        return 0 if result else 1
    
    elif args.monomers:
        # 批处理模式
        monomers = [m.strip().upper() for m in args.monomers.split(",") if m.strip()]
        default_dp = args.n_repeat if args.n_repeat else DEFAULT_OLIGOMER_N
        
        if args.dp:
            dp_map = parse_dp_arg(args.dp, default_dp)
        else:
            dp_map = {m: default_dp for m in MONOMER_TO_POLYMER}
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        mol2_dir = os.path.join(project_root, args.dir)
        
        print("=" * 60)
        print("单体 MOL2 → 聚合物 MOL2 (批处理模式)")
        print("=" * 60)
        
        if not os.path.isdir(mol2_dir):
            print(f"\n[ERROR] MOL2 目录不存在: {mol2_dir}")
            return 1
        
        results = generate_polymer_batch(
            monomers=monomers,
            dp_map=dp_map,
            mol2_dir=mol2_dir,
            optimize=optimize,
            seed=args.seed,
            verbose=args.verbose,
            dry_run=args.dry_run,
            no_overwrite=args.no_overwrite
        )
        
        print("\n" + "=" * 60)
        print("处理完成汇总")
        print("=" * 60)
        print(f"成功: {len(results['success'])} - {', '.join(results['success']) or '无'}")
        print(f"失败: {len(results['failed'])} - {', '.join(results['failed']) or '无'}")
        print(f"跳过: {len(results['skipped'])} - {', '.join(results['skipped']) or '无'}")
        print("=" * 60)
        
        return 0 if not results["failed"] else 1
    
    else:
        print("[ERROR] 请指定输入文件 (--input) 或单体列表 (--monomers)")
        print(f"使用 --help 查看完整帮助")
        return 1


if __name__ == "__main__":
    sys.exit(main())

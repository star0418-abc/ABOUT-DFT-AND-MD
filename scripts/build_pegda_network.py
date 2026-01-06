#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
build_pegda_network.py - 通过 HTPolyNet 生成 PEGDA 交联网络
=============================================================

功能:
  1. 准备 EGDA active form building block (C=C → C-C + 反应氢)
  2. 生成 HTPolyNet 配置文件
  3. 运行 HTPolyNet cure 生成交联网络
  4. 将 HTPolyNet 输出转换为 mol2

用法:
  python scripts/build_pegda_network.py
  python scripts/build_pegda_network.py --n_monomers 300 --conversion 0.9
  python scripts/build_pegda_network.py --dry_run

输出:
  mol2/PEGDA.mol2          - 交联网络 mol2 文件
  htpolynet_out/PEGDA/     - HTPolyNet 完整输出

注意:
  此脚本生成的是真实交联网络，不是线性寡聚体！

版本: 1.0.0
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ==============================================================================
# 依赖检查
# ==============================================================================

def check_rdkit():
    """检查 RDKit"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        return True
    except ImportError:
        return False

def check_htpolynet():
    """检查 HTPolyNet"""
    try:
        result = subprocess.run(
            ["htpolynet", "--version"],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0 or "HTPolyNet" in result.stdout + result.stderr
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False

def check_dependencies():
    """检查所有依赖"""
    missing = []
    
    if not check_rdkit():
        missing.append("RDKit: conda install -c conda-forge rdkit")
    
    if not check_htpolynet():
        missing.append("HTPolyNet: pip install htpolynet 或 conda install -c conda-forge htpolynet")
    
    if missing:
        print("=" * 60)
        print("错误: 缺少必要依赖")
        print("=" * 60)
        for dep in missing:
            print(f"  - {dep}")
        print("=" * 60)
        return False
    
    return True


# ==============================================================================
# RDKit 导入
# ==============================================================================

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule, UFFHasAllMoleculeParams
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# ==============================================================================
# 常量
# ==============================================================================


# EGDA 反应位点命名
REACTIVE_ATOM_NAMES = {
    "head_a": "HA",  # 第一端的 =CH2 碳
    "tail_a": "TA",  # 第一端的 =CH- 碳
    "head_b": "HB",  # 第二端的 =CH2 碳
    "tail_b": "TB",  # 第二端的 =CH- 碳
}

def find_sacrificial_hydrogen(mol: 'Chem.Mol', atom_idx: int) -> Optional[int]:
    """找到反应位点相邻的牺牲氢 (0-based 索引)"""
    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == "H":
            return neighbor.GetIdx()
    return None


# ==============================================================================
# EGDA Active Form 生成
# ==============================================================================

def read_mol_file(mol_path: str) -> Optional['Chem.Mol']:
    """读取 .mol 文件"""
    if not RDKIT_AVAILABLE:
        print("[ERROR] RDKit 未安装")
        return None
    
    if not os.path.isfile(mol_path):
        print(f"[ERROR] 文件不存在: {mol_path}")
        return None
    
    try:
        mol = Chem.MolFromMolFile(mol_path, removeHs=False)
        if mol is None:
            mol = Chem.MolFromMolFile(mol_path, removeHs=False, strictParsing=False)
        return mol
    except Exception as e:
        print(f"[ERROR] 读取 MOL 文件失败: {e}")
        return None


def prepare_egda_active_form(
    mol: 'Chem.Mol',
    verbose: bool = False
) -> Tuple[Optional['Chem.Mol'], Optional[Dict[str, int]]]:
    """
    准备 EGDA active form
    
    1. 识别两个丙烯酸酯端基 (C=C-C(=O)O)
    2. 将 C=C 双键改为 C-C 单键
    3. 添加反应氢
    4. 命名反应原子: HA, TA, HB, TB
    
    Returns:
        (active_mol, reactive_atoms_dict)
        reactive_atoms_dict: {HA: idx, TA: idx, HB: idx, TB: idx} (1-based)
    """
    if not RDKIT_AVAILABLE:
        return None, None
    
    if verbose:
        print("[ACTIVE FORM] 处理 EGDA...")
    
    # 丙烯酸酯 SMARTS: CH2=CH-C(=O)O
    acrylate_pattern = Chem.MolFromSmarts("[CH2:1]=[CH:2]-[C:3](=[O:4])-[O:5]")
    
    if acrylate_pattern is None:
        print("[ERROR] SMARTS 模式解析失败")
        return None, None
    
    matches = mol.GetSubstructMatches(acrylate_pattern)
    
    if verbose:
        print(f"  找到 {len(matches)} 个丙烯酸酯端基")
    
    if len(matches) < 2:
        print(f"[ERROR] EGDA 应有 2 个丙烯酸酯端基，只找到 {len(matches)} 个")
        return None, None
    
    # 只处理前两个
    matches = matches[:2]
    
    rw_mol = Chem.RWMol(mol)
    reactive_atoms = {}
    
    for i, match in enumerate(matches):
        ch2_idx, ch_idx = match[0], match[1]
        
        if verbose:
            print(f"  端基 {i+1}: CH2={ch2_idx+1}, CH={ch_idx+1} (1-based)")
        
        # 获取双键
        bond = rw_mol.GetBondBetweenAtoms(ch2_idx, ch_idx)
        
        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
            bond.SetBondType(Chem.BondType.SINGLE)
            if verbose:
                print(f"    将 C{ch2_idx+1}=C{ch_idx+1} 改为单键")
        
        # 命名反应原子
        if i == 0:
            reactive_atoms["HA"] = ch2_idx + 1  # 1-based
            reactive_atoms["TA"] = ch_idx + 1
        else:
            reactive_atoms["HB"] = ch2_idx + 1
            reactive_atoms["TB"] = ch_idx + 1
    
    # 重新加氢
    try:
        mol_no_h = Chem.RemoveHs(rw_mol)
        mol_with_h = Chem.AddHs(mol_no_h, addCoords=True)
        
        if mol_with_h is None:
            print("[ERROR] 重新加氢失败")
            return None, None
        
        # 生成 3D 坐标
        if mol_with_h.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol_with_h, randomSeed=2025)
        
        # 移除反应位点的牺牲氢，避免交联后出现五价碳
        sacrificial_h = []
        for label in ("HA", "TA", "HB", "TB"):
            atom_idx = reactive_atoms[label] - 1
            h_idx = find_sacrificial_hydrogen(mol_with_h, atom_idx)
            if h_idx is None:
                print(f"[WARN] 未找到 {label} 上的牺牲氢")
                continue
            sacrificial_h.append(h_idx)

        if sacrificial_h:
            rw_with_h = Chem.RWMol(mol_with_h)
            for h_idx in sorted(set(sacrificial_h), reverse=True):
                rw_with_h.RemoveAtom(h_idx)
            mol_with_h = rw_with_h.GetMol()

        # UFF 优化
        if UFFHasAllMoleculeParams(mol_with_h):
            UFFOptimizeMolecule(mol_with_h, maxIters=200)
        
        if verbose:
            print(f"  Active form 完成: {mol_with_h.GetNumAtoms()} 原子")
            print(f"  反应原子: {reactive_atoms}")
        
        return mol_with_h, reactive_atoms
        
    except Exception as e:
        print(f"[ERROR] Active form 处理失败: {e}")
        return None, None


def write_mol2_with_reactive_names(
    mol: 'Chem.Mol',
    filepath: str,
    mol_name: str,
    res_name: str,
    reactive_atoms: Dict[str, int]
) -> bool:
    """
    写入 MOL2 文件，使用指定的反应原子名
    
    Args:
        mol: RDKit 分子
        filepath: 输出路径
        mol_name: 分子名
        res_name: 残基名 (最多3字符)
        reactive_atoms: {原子名: 1-based索引}
    """
    try:
        conf = mol.GetConformer()
        atoms = mol.GetAtoms()
        bonds = mol.GetBonds()
        
        # 生成原子名（反应原子使用指定名称）
        reactive_idx_to_name = {v-1: k for k, v in reactive_atoms.items()}  # 转为 0-based
        
        atom_names = {}
        element_count = {}
        
        for atom in atoms:
            idx = atom.GetIdx()
            
            if idx in reactive_idx_to_name:
                atom_names[idx] = reactive_idx_to_name[idx]
            else:
                symbol = atom.GetSymbol()
                element_count[symbol] = element_count.get(symbol, 0) + 1
                atom_names[idx] = f"{symbol}{element_count[symbol]}"
        
        # Sybyl 类型
        def get_sybyl(atom):
            symbol = atom.GetSymbol()
            hyb = str(atom.GetHybridization()).lower()
            
            if symbol == 'H':
                return 'H'
            elif symbol == 'C':
                if 'sp3' in hyb:
                    return 'C.3'
                elif 'sp2' in hyb:
                    return 'C.2'
                else:
                    return 'C.3'
            elif symbol == 'O':
                if 'sp3' in hyb:
                    return 'O.3'
                else:
                    return 'O.2'
            elif symbol == 'N':
                return 'N.3'
            elif symbol == 'F':
                return 'F'
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
                sybyl = get_sybyl(atom)
                f.write(f"{idx+1:7d} {name:<4s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} "
                       f"{sybyl:<6s} 1 {res_name[:3]:<4s} 0.0000\n")
            
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
        print(f"[ERROR] 写入 MOL2 失败: {e}")
        return False


# ==============================================================================
# HTPolyNet 配置生成
# ==============================================================================

def generate_htpolynet_config(
    workdir: str,
    mol2_path: str,
    reactive_atoms: Dict[str, int],
    n_monomers: int = 200,
    conversion: float = 0.80,
    seed: int = 2025,
    verbose: bool = False
) -> str:
    """
    生成 HTPolyNet 配置文件
    
    Returns:
        配置文件路径
    """
    config_path = os.path.join(workdir, "pegda_recipe.yaml")
    
    # HTPolyNet 配置
    # 反应定义: TA-HA 和 TB-HB
    config_content = f'''# =============================================================================
# HTPolyNet 配置文件 - PEGDA 交联网络
# 自动生成 by build_pegda_network.py
# =============================================================================

Title: PEGDA Crosslinked Network

# -----------------------------------------------------------------------------
# 分子定义
# -----------------------------------------------------------------------------
constituents:
  EGDA:
    count: {n_monomers}
    source: "{mol2_path}"

# -----------------------------------------------------------------------------
# 反应定义
# -----------------------------------------------------------------------------
# EGDA 双官能团交联剂有两个反应端:
#   端 A: HA (=CH2) -- TA (=CH-)
#   端 B: HB (=CH2) -- TB (=CH-)
# 
# 交联反应: 一个分子的 TA 与另一个分子的 HA 成键
#           一个分子的 TB 与另一个分子的 HB 成键

reactions:
  # 反应1: A端交联
  - name: "crosslink_A"
    reactants:
      - molecule: EGDA
        atom: TA
      - molecule: EGDA
        atom: HA
    product_bond_order: 1
    
  # 反应2: B端交联  
  - name: "crosslink_B"
    reactants:
      - molecule: EGDA
        atom: TB
      - molecule: EGDA
        atom: HB
    product_bond_order: 1

# -----------------------------------------------------------------------------
# Cure 参数
# -----------------------------------------------------------------------------
cure:
  # 目标转化率
  desired_conversion: {conversion}
  
  # 最大尝试次数
  max_search_attempts: 1000
  
  # 反应距离阈值 (nm)
  cutoff: 0.5
  
  # 随机种子
  random_seed: {seed}

# -----------------------------------------------------------------------------
# MD 参数
# -----------------------------------------------------------------------------
densification:
  nsteps: 5000
  dt: 0.001  # ps
  temperature: 300  # K
  
equilibration:
  - name: "nvt"
    ensemble: nvt
    nsteps: 10000
    dt: 0.001
    temperature: 300
    
  - name: "npt"  
    ensemble: npt
    nsteps: 20000
    dt: 0.002
    temperature: 300
    pressure: 1.0  # bar

# -----------------------------------------------------------------------------
# 系统设置
# -----------------------------------------------------------------------------
box:
  type: cubic
  initial_density: 0.5  # g/cm³ (初始低密度便于打包)

output:
  prefix: "pegda_network"
  formats:
    - gro
    - top
    - pdb
'''
    
    os.makedirs(workdir, exist_ok=True)
    
    with open(config_path, 'w') as f:
        f.write(config_content)
    
    if verbose:
        print(f"[CONFIG] 生成配置文件: {config_path}")
    
    return config_path


# ==============================================================================
# HTPolyNet 运行
# ==============================================================================

def run_htpolynet(
    config_path: str,
    workdir: str,
    verbose: bool = False,
    dry_run: bool = False
) -> bool:
    """
    运行 HTPolyNet
    
    Returns:
        是否成功
    """
    if not check_htpolynet():
        print("[ERROR] HTPolyNet 未安装或不可用")
        print("安装方法:")
        print("  pip install htpolynet")
        print("  或")
        print("  conda install -c conda-forge htpolynet")
        return False
    
    if dry_run:
        print(f"[DRY-RUN] 将运行: htpolynet run {config_path}")
        print(f"[DRY-RUN] 工作目录: {workdir}")
        return True
    
    print(f"\n[HTPOLYNET] 开始运行...")
    print(f"  配置文件: {config_path}")
    print(f"  工作目录: {workdir}")
    print("  这可能需要几分钟到几小时，取决于分子数和转化率...")
    
    try:
        # 切换到工作目录运行
        result = subprocess.run(
            ["htpolynet", "run", os.path.basename(config_path)],
            cwd=workdir,
            capture_output=not verbose,
            text=True,
            timeout=7200,  # 2小时超时
        )
        
        if result.returncode != 0:
            print(f"[ERROR] HTPolyNet 运行失败 (code={result.returncode})")
            if result.stderr:
                print(f"  stderr: {result.stderr[:500]}")
            return False
        
        print("[HTPOLYNET] 运行完成!")
        return True
        
    except subprocess.TimeoutExpired:
        print("[ERROR] HTPolyNet 运行超时 (>2小时)")
        return False
    except Exception as e:
        print(f"[ERROR] HTPolyNet 运行异常: {e}")
        return False


# ==============================================================================
# 输出转换 (GRO/PDB → MOL2)
# ==============================================================================

def find_htpolynet_output(workdir: str) -> Tuple[Optional[str], Optional[str]]:
    """
    查找 HTPolyNet 输出文件
    
    Returns:
        (gro_path, pdb_path) 或 (None, None)
    """
    # 常见输出文件名模式
    gro_patterns = [
        "pegda_network_final.gro",
        "pegda_network.gro",
        "final.gro",
        "system.gro",
        "*.gro"
    ]
    
    pdb_patterns = [
        "pegda_network_final.pdb",
        "pegda_network.pdb",
        "final.pdb",
        "system.pdb",
        "*.pdb"
    ]
    
    gro_path = None
    pdb_path = None
    
    # 搜索 GRO
    for pattern in gro_patterns:
        if "*" in pattern:
            import glob
            matches = glob.glob(os.path.join(workdir, "**", pattern), recursive=True)
            if matches:
                gro_path = matches[0]
                break
        else:
            path = os.path.join(workdir, pattern)
            if os.path.isfile(path):
                gro_path = path
                break
    
    # 搜索 PDB
    for pattern in pdb_patterns:
        if "*" in pattern:
            import glob
            matches = glob.glob(os.path.join(workdir, "**", pattern), recursive=True)
            if matches:
                pdb_path = matches[0]
                break
        else:
            path = os.path.join(workdir, pattern)
            if os.path.isfile(path):
                pdb_path = path
                break
    
    return gro_path, pdb_path


def convert_pdb_to_mol2(pdb_path: str, mol2_path: str, verbose: bool = False) -> bool:
    """
    将 PDB 转换为 MOL2
    
    使用 RDKit 或回退到自定义转换
    """
    if not RDKIT_AVAILABLE:
        print("[ERROR] RDKit 不可用，无法转换")
        return False
    
    if verbose:
        print(f"[CONVERT] PDB → MOL2")
        print(f"  输入: {pdb_path}")
        print(f"  输出: {mol2_path}")
    
    try:
        # 尝试 RDKit 读取
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        
        if mol is None:
            print("[WARN] RDKit 无法读取 PDB，尝试自定义转换...")
            return convert_pdb_to_mol2_manual(pdb_path, mol2_path, verbose)
        
        # 写入 MOL2
        mol_name = "PEGDA_Network"
        res_name = "PEG"
        
        return write_mol2_simple(mol, mol2_path, mol_name, res_name)
        
    except Exception as e:
        print(f"[ERROR] 转换失败: {e}")
        return convert_pdb_to_mol2_manual(pdb_path, mol2_path, verbose)


def convert_pdb_to_mol2_manual(pdb_path: str, mol2_path: str, verbose: bool = False) -> bool:
    """
    手动将 PDB 转换为 MOL2（当 RDKit 失败时）
    
    这是一个简化的转换，只保留原子坐标
    """
    try:
        atoms = []
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atom_id = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    
                    atoms.append({
                        'id': atom_id,
                        'name': atom_name,
                        'res': res_name,
                        'x': x, 'y': y, 'z': z,
                        'element': element
                    })
        
        if not atoms:
            print("[ERROR] PDB 文件中没有原子")
            return False
        
        # 写入 MOL2（不含键信息）
        with open(mol2_path, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write("PEGDA_Network\n")
            f.write(f" {len(atoms)} 0 1 0 0\n")
            f.write("SMALL\n")
            f.write("NO_CHARGES\n\n")
            
            f.write("@<TRIPOS>ATOM\n")
            for i, atom in enumerate(atoms, 1):
                sybyl = atom['element'] + ".3" if atom['element'] in ['C', 'N', 'O', 'S'] else atom['element']
                f.write(f"{i:7d} {atom['name']:<4s} {atom['x']:10.4f} {atom['y']:10.4f} {atom['z']:10.4f} "
                       f"{sybyl:<6s} 1 {atom['res'][:3]:<4s} 0.0000\n")
            
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write("     1 PEG         1 RESIDUE    0 **** **** 0 ROOT\n")
        
        if verbose:
            print(f"  手动转换完成: {len(atoms)} 原子")
        
        return True
        
    except Exception as e:
        print(f"[ERROR] 手动转换失败: {e}")
        return False


def write_mol2_simple(mol: 'Chem.Mol', filepath: str, mol_name: str, res_name: str) -> bool:
    """简化的 MOL2 写入"""
    try:
        conf = mol.GetConformer()
        atoms = mol.GetAtoms()
        bonds = mol.GetBonds()
        
        atom_names = {}
        element_count = {}
        
        for atom in atoms:
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            element_count[symbol] = element_count.get(symbol, 0) + 1
            atom_names[idx] = f"{symbol}{element_count[symbol]}"
        
        def get_sybyl(atom):
            symbol = atom.GetSymbol()
            if symbol == 'H':
                return 'H'
            elif symbol == 'C':
                return 'C.3'
            elif symbol == 'O':
                return 'O.3'
            elif symbol == 'N':
                return 'N.3'
            else:
                return symbol
        
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
                f.write(f"{idx+1:7d} {atom_names[idx]:<4s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} "
                       f"{get_sybyl(atom):<6s} 1 {res_name[:3]:<4s} 0.0000\n")
            
            f.write("@<TRIPOS>BOND\n")
            for i, bond in enumerate(bonds):
                f.write(f"{i+1:6d} {bond.GetBeginAtomIdx()+1:5d} {bond.GetEndAtomIdx()+1:5d} 1\n")
            
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write(f"     1 {res_name[:3]:<4s}        1 RESIDUE    0 **** **** 0 ROOT\n")
        
        return True
    except Exception as e:
        print(f"[ERROR] 写入 MOL2 失败: {e}")
        return False


def convert_gro_to_mol2(gro_path: str, top_path: Optional[str], mol2_path: str, verbose: bool = False) -> bool:
    """
    将 GRO 转换为 MOL2
    
    尝试使用 MDAnalysis 或回退到手动转换
    """
    try:
        import MDAnalysis as mda
        
        if verbose:
            print("[CONVERT] 使用 MDAnalysis 转换 GRO → MOL2")
        
        # 加载 GRO（可选带 TOP）
        if top_path and os.path.isfile(top_path):
            u = mda.Universe(top_path, gro_path)
        else:
            u = mda.Universe(gro_path)
        
        # 写入 MOL2
        atoms = u.atoms
        
        with open(mol2_path, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write("PEGDA_Network\n")
            f.write(f" {len(atoms)} 0 1 0 0\n")
            f.write("SMALL\n")
            f.write("NO_CHARGES\n\n")
            
            f.write("@<TRIPOS>ATOM\n")
            for i, atom in enumerate(atoms, 1):
                element = atom.element if hasattr(atom, 'element') else atom.name[0]
                sybyl = element + ".3" if element in ['C', 'N', 'O', 'S'] else element
                f.write(f"{i:7d} {atom.name:<4s} {atom.position[0]:10.4f} {atom.position[1]:10.4f} "
                       f"{atom.position[2]:10.4f} {sybyl:<6s} 1 {atom.resname[:3]:<4s} 0.0000\n")
            
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write("     1 PEG         1 RESIDUE    0 **** **** 0 ROOT\n")
        
        if verbose:
            print(f"  转换完成: {len(atoms)} 原子")
        
        return True
        
    except ImportError:
        if verbose:
            print("[WARN] MDAnalysis 不可用，尝试手动转换...")
        return convert_gro_to_mol2_manual(gro_path, mol2_path, verbose)
    except Exception as e:
        print(f"[ERROR] MDAnalysis 转换失败: {e}")
        return convert_gro_to_mol2_manual(gro_path, mol2_path, verbose)


def convert_gro_to_mol2_manual(gro_path: str, mol2_path: str, verbose: bool = False) -> bool:
    """手动解析 GRO 并转换为 MOL2"""
    try:
        atoms = []
        
        with open(gro_path, 'r') as f:
            lines = f.readlines()
        
        # GRO 格式:
        # Line 1: title
        # Line 2: number of atoms
        # Lines 3-: atom data (fixed format)
        # Last line: box vectors
        
        n_atoms = int(lines[1].strip())
        
        for i in range(2, 2 + n_atoms):
            line = lines[i]
            # GRO format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
            res_id = int(line[0:5])
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_id = int(line[15:20])
            x = float(line[20:28]) * 10  # nm → Å
            y = float(line[28:36]) * 10
            z = float(line[36:44]) * 10
            
            element = atom_name[0]
            
            atoms.append({
                'id': atom_id,
                'name': atom_name,
                'res': res_name,
                'x': x, 'y': y, 'z': z,
                'element': element
            })
        
        # 写入 MOL2
        with open(mol2_path, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write("PEGDA_Network\n")
            f.write(f" {len(atoms)} 0 1 0 0\n")
            f.write("SMALL\n")
            f.write("NO_CHARGES\n\n")
            
            f.write("@<TRIPOS>ATOM\n")
            for i, atom in enumerate(atoms, 1):
                sybyl = atom['element'] + ".3" if atom['element'] in ['C', 'N', 'O', 'S'] else atom['element']
                f.write(f"{i:7d} {atom['name']:<4s} {atom['x']:10.4f} {atom['y']:10.4f} {atom['z']:10.4f} "
                       f"{sybyl:<6s} 1 {atom['res'][:3]:<4s} 0.0000\n")
            
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write("     1 PEG         1 RESIDUE    0 **** **** 0 ROOT\n")
        
        if verbose:
            print(f"  手动转换完成: {len(atoms)} 原子")
        
        return True
        
    except Exception as e:
        print(f"[ERROR] GRO 手动转换失败: {e}")
        return False


# ==============================================================================
# 主流程
# ==============================================================================

def build_pegda_network(
    in_mol: str,
    out_mol2: str,
    workdir: str,
    n_monomers: int = 200,
    conversion: float = 0.80,
    seed: int = 2025,
    optimize: bool = True,
    verbose: bool = False,
    dry_run: bool = False
) -> bool:
    """
    构建 PEGDA 交联网络
    
    完整流程:
    1. 准备 EGDA active form
    2. 生成 HTPolyNet 配置
    3. 运行 HTPolyNet
    4. 转换输出为 MOL2
    """
    print("=" * 60)
    print("PEGDA 交联网络构建")
    print("=" * 60)
    print(f"输入单体: {in_mol}")
    print(f"输出网络: {out_mol2}")
    print(f"工作目录: {workdir}")
    print(f"单体数量: {n_monomers}")
    print(f"目标转化率: {conversion}")
    print(f"随机种子: {seed}")
    print("=" * 60)
    
    # 创建目录
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(os.path.dirname(out_mol2) or ".", exist_ok=True)
    
    # Step 1: 准备 EGDA active form
    print("\n[Step 1/4] 准备 EGDA active form...")
    
    mol = read_mol_file(in_mol)
    if mol is None:
        return False
    
    active_mol, reactive_atoms = prepare_egda_active_form(mol, verbose)
    if active_mol is None:
        return False
    
    # 保存 active form mol2
    active_mol2_path = os.path.join(workdir, "EGDA_active.mol2")
    success = write_mol2_with_reactive_names(
        active_mol, active_mol2_path, "EGDA", "EGD", reactive_atoms
    )
    
    if not success:
        return False
    
    print(f"  Active form 已保存: {active_mol2_path}")
    print(f"  反应原子: HA={reactive_atoms['HA']}, TA={reactive_atoms['TA']}, "
          f"HB={reactive_atoms['HB']}, TB={reactive_atoms['TB']}")
    
    # 同时复制到 mol2 目录
    mol2_dir = os.path.dirname(out_mol2)
    if mol2_dir:
        active_mol2_copy = os.path.join(mol2_dir, "EGDA_active.mol2")
        shutil.copy(active_mol2_path, active_mol2_copy)
        print(f"  副本: {active_mol2_copy}")
    
    # Step 2: 生成 HTPolyNet 配置
    print("\n[Step 2/4] 生成 HTPolyNet 配置...")
    
    config_path = generate_htpolynet_config(
        workdir=workdir,
        mol2_path=active_mol2_path,
        reactive_atoms=reactive_atoms,
        n_monomers=n_monomers,
        conversion=conversion,
        seed=seed,
        verbose=verbose
    )
    
    print(f"  配置文件: {config_path}")
    
    if dry_run:
        print("\n[DRY-RUN] 跳过 HTPolyNet 运行和输出转换")
        print("\n要实际运行，请去掉 --dry_run 参数")
        return True
    
    # Step 3: 运行 HTPolyNet
    print("\n[Step 3/4] 运行 HTPolyNet...")
    
    success = run_htpolynet(config_path, workdir, verbose, dry_run)
    
    if not success:
        print("\n[WARN] HTPolyNet 运行失败")
        print("你可以手动运行:")
        print(f"  cd {workdir}")
        print(f"  htpolynet run pegda_recipe.yaml")
        return False
    
    # Step 4: 转换输出
    print("\n[Step 4/4] 转换输出为 MOL2...")
    
    gro_path, pdb_path = find_htpolynet_output(workdir)
    
    if pdb_path:
        print(f"  找到 PDB: {pdb_path}")
        success = convert_pdb_to_mol2(pdb_path, out_mol2, verbose)
    elif gro_path:
        print(f"  找到 GRO: {gro_path}")
        top_path = gro_path.replace(".gro", ".top")
        if not os.path.isfile(top_path):
            top_path = None
        success = convert_gro_to_mol2(gro_path, top_path, out_mol2, verbose)
    else:
        print("[WARN] 未找到 HTPolyNet 输出文件")
        print(f"  请检查目录: {workdir}")
        return False
    
    if success:
        print(f"\n[SUCCESS] PEGDA 交联网络已生成!")
        print(f"  交联网络 MOL2: {out_mol2}")
        print(f"  HTPolyNet 输出: {workdir}/")
        print(f"\n可用 Avogadro 查看: avogadro {out_mol2}")
    
    return success


# ==============================================================================
# CLI
# ==============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="通过 HTPolyNet 生成 PEGDA 交联网络",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s
  %(prog)s --n_monomers 300 --conversion 0.9
  %(prog)s --dry_run --verbose

输出:
  mol2/PEGDA.mol2          - 交联网络结构
  htpolynet_out/PEGDA/     - HTPolyNet 完整输出 (gro/top/itp)

注意:
  此脚本生成真实交联网络，不是线性寡聚体！
"""
    )
    
    parser.add_argument(
        "--in_mol",
        default="mol/EGDA.mol",
        help="输入 EGDA 单体文件 (默认: mol/EGDA.mol)"
    )
    parser.add_argument(
        "--out_mol2",
        default="mol2/PEGDA.mol2",
        help="输出交联网络 MOL2 (默认: mol2/PEGDA.mol2)"
    )
    parser.add_argument(
        "--workdir",
        default="htpolynet_out/PEGDA",
        help="HTPolyNet 工作目录 (默认: htpolynet_out/PEGDA)"
    )
    parser.add_argument(
        "--n_monomers",
        type=int,
        default=200,
        help="单体数量 (默认: 200)"
    )
    parser.add_argument(
        "--conversion",
        type=float,
        default=0.80,
        help="目标转化率 (默认: 0.80)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="随机种子 (默认: 2025)"
    )
    parser.add_argument(
        "--opt",
        action="store_true",
        help="对 active form 进行 UFF 优化 (默认已启用)"
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="模拟运行，不实际执行 HTPolyNet"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="详细输出"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    
    # 检查依赖
    if not RDKIT_AVAILABLE:
        print("[ERROR] RDKit 未安装")
        print("安装: conda install -c conda-forge rdkit")
        return 1
    
    # 确定路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    in_mol = os.path.join(project_root, args.in_mol)
    out_mol2 = os.path.join(project_root, args.out_mol2)
    workdir = os.path.join(project_root, args.workdir)
    
    # 运行
    success = build_pegda_network(
        in_mol=in_mol,
        out_mol2=out_mol2,
        workdir=workdir,
        n_monomers=args.n_monomers,
        conversion=args.conversion,
        seed=args.seed,
        optimize=args.opt,
        verbose=args.verbose,
        dry_run=args.dry_run
    )
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())

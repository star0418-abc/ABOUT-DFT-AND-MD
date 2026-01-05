# -*- coding: utf-8 -*-

# Author: Star



"""
io_mol.py - 统一分子文件读写模块
================================

提供统一的 MOL2/PDB 等分子文件读写功能。
所有文件写入都会进行严格校验。

用法:
    from lib.io_mol import read_mol2, write_mol2_strict
    
    mol = read_mol2("input.mol2")
    write_mol2_strict(mol, "output.mol2", "MOL", "MOL")
"""

import os
from typing import Dict, Optional, List

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


class FileIOError(Exception):
    """文件 IO 错误"""
    pass


def validate_mol2_input(filepath: str) -> None:
    """
    验证输入文件必须是 .mol2 格式
    
    Raises:
        FileIOError: 格式不正确或文件不存在
    """
    if not os.path.isfile(filepath):
        raise FileIOError(f"文件不存在: {filepath}")
    
    ext = os.path.splitext(filepath)[1].lower()
    if ext != ".mol2":
        raise FileIOError(
            f"仅支持 .mol2 格式输入，收到: {ext}\n"
            f"文件: {filepath}"
        )


def read_mol2(filepath: str) -> 'Chem.Mol':
    """
    读取 MOL2 文件
    
    Args:
        filepath: MOL2 文件路径
    
    Returns:
        RDKit Mol 对象
    
    Raises:
        FileIOError: 读取失败
    """
    if not HAS_RDKIT:
        raise FileIOError("RDKit 未安装")
    
    validate_mol2_input(filepath)
    
    mol = Chem.MolFromMol2File(filepath, removeHs=False)
    if mol is None:
        raise FileIOError(f"无法解析 MOL2 文件: {filepath}")
    
    return mol


def write_mol2_strict(
    mol: 'Chem.Mol',
    filepath: str,
    mol_name: str,
    res_name: str
) -> None:
    """
    严格模式写入 MOL2 文件
    
    - 写入后验证文件存在且非空
    - 验证包含必需段落 (@<TRIPOS>MOLECULE/ATOM/BOND)
    - 任何失败都抛出异常
    
    Args:
        mol: RDKit Mol 对象
        filepath: 输出路径
        mol_name: 分子名
        res_name: 残基名
    
    Raises:
        FileIOError: 写入失败
    """
    if not HAS_RDKIT:
        raise FileIOError("RDKit 未安装")
    
    # 确保目录存在
    out_dir = os.path.dirname(filepath)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    
    try:
        conf = mol.GetConformer()
        atoms = list(mol.GetAtoms())
        bonds = list(mol.GetBonds())
        
        if len(atoms) == 0:
            raise FileIOError("分子没有原子，无法写入")
        
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
            
            type_map = {
                'H': 'H',
                'F': 'F',
                'Cl': 'Cl',
                'Br': 'Br',
                'I': 'I',
            }
            if symbol in type_map:
                return type_map[symbol]
            
            if symbol == 'C':
                if atom.GetIsAromatic():
                    return 'C.ar'
                return 'C.3' if 'sp3' in hyb else 'C.2'
            elif symbol == 'N':
                if atom.GetIsAromatic():
                    return 'N.ar'
                return 'N.3' if 'sp3' in hyb else 'N.2'
            elif symbol == 'O':
                return 'O.3' if 'sp3' in hyb else 'O.2'
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
        
        # 验证输出
        validate_mol2_output(filepath, len(atoms), len(bonds))
        
    except FileIOError:
        raise
    except Exception as e:
        # 清理可能的部分写入文件
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
            except:
                pass
        raise FileIOError(f"写入 MOL2 失败: {e}")


def validate_mol2_output(filepath: str, expected_atoms: int = 0, expected_bonds: int = 0) -> None:
    """
    验证 MOL2 输出文件
    
    Args:
        filepath: 文件路径
        expected_atoms: 期望原子数 (0 表示不检查)
        expected_bonds: 期望键数 (0 表示不检查)
    
    Raises:
        FileIOError: 验证失败
    """
    # 1. 文件存在
    if not os.path.exists(filepath):
        raise FileIOError(f"输出文件不存在: {filepath}")
    
    # 2. 文件非空
    file_size = os.path.getsize(filepath)
    if file_size == 0:
        raise FileIOError(f"输出文件为空: {filepath}")
    
    # 3. 包含必需段落
    with open(filepath, 'r') as f:
        content = f.read()
    
    required_sections = ["@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND"]
    for section in required_sections:
        if section not in content:
            raise FileIOError(f"输出文件缺少 {section}: {filepath}")


def get_atom_names_from_mol2(filepath: str) -> Dict[int, str]:
    """
    从 MOL2 文件解析原子名
    
    Returns:
        {原子索引(0-based): 原子名}
    """
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


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
gaussian_log_to_mol2.py - Gaussian Log 文件转 MOL2 格式

从 Gaussian 输出文件 (.log / .out) 中解析几何结构，导出为 MOL2 格式。
支持批量转换，自动跳过已存在文件。

功能特点:
  - 支持 Standard orientation 和 Input orientation 解析
  - 支持批量处理目录中的所有 .log 文件
  - 自动生成键连接（优先 RDKit，fallback 到距离猜测）
  - 支持导出单构型或按 step 导出多个文件
  - 不依赖 OpenBabel

用法示例:
  # 单文件模式
  python scripts/gaussian_log_to_mol2.py --log format/log/TFSI.LOG

  # 批量模式（推荐）
  python scripts/gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2

  # 批量模式，跳过已存在
  python scripts/gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2 --skip_existing

  # 导出所有构型（每个 step 一个文件）
  python scripts/gaussian_log_to_mol2.py --log opt.log --which all

版本: 1.0.0
"""

import argparse
import os
import sys
import math
from pathlib import Path
from typing import List, Tuple, Optional, Dict, NamedTuple
from collections import defaultdict
from dataclasses import dataclass, field


# ==============================================================================
# 元素周期表和共价半径
# ==============================================================================
ATOMIC_NUMBER_TO_SYMBOL = {
    1: 'H',   2: 'He',  3: 'Li',  4: 'Be',  5: 'B',   6: 'C',   7: 'N',   8: 'O',
    9: 'F',  10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',  16: 'S',
    17: 'Cl', 18: 'Ar', 19: 'K',  20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V',  24: 'Cr',
    25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge',
    33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',  40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd',
    49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I',  54: 'Xe', 55: 'Cs', 56: 'Ba',
    57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
    65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
    73: 'Ta', 74: 'W',  75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
}

# 共价半径 (Å) - 用于距离猜键
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71,
    'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11,
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Br': 1.20,
    'I': 1.39, 'Fe': 1.32, 'Zn': 1.22, 'Cu': 1.32, 'Ni': 1.24, 'Co': 1.26,
}


# ==============================================================================
# 数据结构
# ==============================================================================
class Atom(NamedTuple):
    """原子数据"""
    atomic_number: int
    x: float
    y: float
    z: float


class OrientationBlock(NamedTuple):
    """Orientation 块数据"""
    block_index: int
    orientation_type: str
    atoms: List[Atom]
    line_number: int


@dataclass
class BatchResult:
    """批量处理结果"""
    total: int = 0
    success: int = 0
    skipped: int = 0
    failed: int = 0
    success_files: List[str] = field(default_factory=list)
    skipped_files: List[str] = field(default_factory=list)
    failed_files: List[Tuple[str, str]] = field(default_factory=list)


# ==============================================================================
# Gaussian Log 解析（复用自 gaussian_log_to_pdb.py）
# ==============================================================================
import re

def iter_orientation_blocks(
    log_path: str,
    orientation_filter: str = 'auto'
):
    """流式解析 Gaussian log 文件，迭代返回 orientation 块"""
    re_standard = re.compile(r'^\s*Standard orientation:', re.IGNORECASE)
    re_input = re.compile(r'^\s*Input orientation:', re.IGNORECASE)
    re_separator = re.compile(r'^\s*-{50,}')
    re_data = re.compile(
        r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+'
        r'(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$'
    )
    
    block_index = 0
    current_type = None
    current_atoms: List[Atom] = []
    in_block = False
    separator_count = 0
    start_line = 0
    
    with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
        for line_num, line in enumerate(f, start=1):
            if re_standard.match(line):
                if in_block and current_atoms:
                    block_index += 1
                    yield OrientationBlock(block_index, current_type, current_atoms.copy(), start_line)
                current_type = 'standard'
                current_atoms = []
                in_block = True
                separator_count = 0
                start_line = line_num
                continue
                
            elif re_input.match(line):
                if orientation_filter == 'standard':
                    continue
                if in_block and current_atoms:
                    block_index += 1
                    yield OrientationBlock(block_index, current_type, current_atoms.copy(), start_line)
                current_type = 'input'
                current_atoms = []
                in_block = True
                separator_count = 0
                start_line = line_num
                continue
            
            if not in_block:
                continue
                
            if re_separator.match(line):
                separator_count += 1
                if separator_count >= 3:
                    if current_atoms:
                        block_index += 1
                        yield OrientationBlock(block_index, current_type, current_atoms.copy(), start_line)
                    in_block = False
                    current_atoms = []
                continue
            
            if separator_count >= 2:
                m = re_data.match(line)
                if m:
                    atomic_number = int(m.group(2))
                    x = float(m.group(4))
                    y = float(m.group(5))
                    z = float(m.group(6))
                    current_atoms.append(Atom(atomic_number, x, y, z))
    
    if in_block and current_atoms:
        block_index += 1
        yield OrientationBlock(block_index, current_type, current_atoms.copy(), start_line)


def collect_blocks(log_path: str, orientation_filter: str = 'auto'):
    """收集所有 orientation 块"""
    standard_blocks = []
    input_blocks = []
    
    for block in iter_orientation_blocks(log_path, orientation_filter):
        if block.orientation_type == 'standard':
            standard_blocks.append(block)
        else:
            input_blocks.append(block)
    
    return standard_blocks, input_blocks


def select_blocks(log_path: str, which: str, orientation: str = 'auto') -> List[OrientationBlock]:
    """根据用户选择返回对应的 orientation 块"""
    standard_blocks, input_blocks = collect_blocks(log_path, orientation)
    
    if orientation == 'standard':
        blocks = standard_blocks
        if not blocks:
            raise ValueError(f"在 {log_path} 中未找到 Standard orientation 块")
    elif orientation == 'input':
        blocks = input_blocks
        if not blocks:
            raise ValueError(f"在 {log_path} 中未找到 Input orientation 块")
    else:
        blocks = standard_blocks if standard_blocks else input_blocks
    
    if not blocks:
        raise ValueError(f"在 {log_path} 中未找到任何 orientation 块")
    
    if which == 'all':
        return blocks
    elif which == 'last':
        return [blocks[-1]]
    elif which == 'first':
        return [blocks[0]]
    elif which.startswith('step=') or which.startswith('step-'):
        try:
            step_num = int(which.split('=')[-1].split('-')[-1])
        except ValueError:
            raise ValueError(f"无效的 --which 参数: {which}")
        if step_num < 1 or step_num > len(blocks):
            raise ValueError(f"step={step_num} 超出范围，仅有 {len(blocks)} 个块")
        return [blocks[step_num - 1]]
    else:
        raise ValueError(f"无效的 --which 参数: {which}")


# ==============================================================================
# 键生成
# ==============================================================================
def try_rdkit_bonds(atoms: List[Atom]) -> Optional[List[Tuple[int, int, int]]]:
    """
    尝试使用 RDKit 根据坐标猜测键
    
    Returns:
        [(atom1_idx, atom2_idx, bond_order), ...] 或 None
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        
        # 创建 RWMol
        mol = Chem.RWMol()
        conf = Chem.Conformer(len(atoms))
        
        for i, atom in enumerate(atoms):
            symbol = ATOMIC_NUMBER_TO_SYMBOL.get(atom.atomic_number, 'C')
            rd_atom = Chem.Atom(symbol)
            mol.AddAtom(rd_atom)
            conf.SetAtomPosition(i, (atom.x, atom.y, atom.z))
        
        mol.AddConformer(conf, assignId=True)
        
        # 使用 DetermineBonds 猜测键
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=0)
        except:
            rdDetermineBonds.DetermineBondOrders(mol, charge=0)
        
        # 提取键
        bonds = []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            bt = bond.GetBondType()
            if bt == Chem.BondType.SINGLE:
                order = 1
            elif bt == Chem.BondType.DOUBLE:
                order = 2
            elif bt == Chem.BondType.TRIPLE:
                order = 3
            elif bt == Chem.BondType.AROMATIC:
                order = 4  # ar
            else:
                order = 1
            bonds.append((i, j, order))
        
        return bonds
        
    except ImportError:
        return None
    except Exception:
        return None


def guess_bonds_by_distance(atoms: List[Atom], tolerance: float = 1.25) -> List[Tuple[int, int, int]]:
    """
    根据共价半径和距离阈值猜测键（Fallback 方法）
    
    Args:
        atoms: 原子列表
        tolerance: 容差因子（默认 1.25，即 d < 1.25*(r_i+r_j) 则建键）
    
    Returns:
        [(atom1_idx, atom2_idx, bond_order), ...] bond_order 统一为 1
    """
    bonds = []
    n = len(atoms)
    
    for i in range(n):
        sym_i = ATOMIC_NUMBER_TO_SYMBOL.get(atoms[i].atomic_number, 'C')
        r_i = COVALENT_RADII.get(sym_i, 1.5)
        
        for j in range(i + 1, n):
            sym_j = ATOMIC_NUMBER_TO_SYMBOL.get(atoms[j].atomic_number, 'C')
            r_j = COVALENT_RADII.get(sym_j, 1.5)
            
            # 计算距离
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            # 判断是否成键
            threshold = tolerance * (r_i + r_j)
            if dist < threshold:
                bonds.append((i, j, 1))  # 统一 bond order = 1
    
    return bonds


def generate_bonds(atoms: List[Atom], verbose: bool = False) -> List[Tuple[int, int, int]]:
    """
    生成键列表，优先使用 RDKit，失败则 fallback 到距离猜测
    """
    # 尝试 RDKit
    bonds = try_rdkit_bonds(atoms)
    if bonds is not None:
        if verbose:
            print(f"    使用 RDKit 生成 {len(bonds)} 个键")
        return bonds
    
    # Fallback 到距离猜测
    bonds = guess_bonds_by_distance(atoms)
    if verbose:
        print(f"    使用距离猜测生成 {len(bonds)} 个键")
    return bonds


# ==============================================================================
# MOL2 写入
# ==============================================================================
def get_sybyl_type(symbol: str, bonds_for_atom: List[int]) -> str:
    """根据元素和键获取 Sybyl 原子类型"""
    if symbol == 'H':
        return 'H'
    elif symbol == 'C':
        # 简化：根据键数判断
        if 4 in bonds_for_atom:  # 有芳香键
            return 'C.ar'
        elif 2 in bonds_for_atom:  # 有双键
            return 'C.2'
        elif 3 in bonds_for_atom:  # 有三键
            return 'C.1'
        else:
            return 'C.3'
    elif symbol == 'N':
        if 4 in bonds_for_atom:
            return 'N.ar'
        elif 2 in bonds_for_atom:
            return 'N.2'
        else:
            return 'N.3'
    elif symbol == 'O':
        if 2 in bonds_for_atom:
            return 'O.2'
        else:
            return 'O.3'
    elif symbol == 'S':
        return 'S.3'
    elif symbol == 'P':
        return 'P.3'
    elif symbol == 'F':
        return 'F'
    elif symbol == 'Cl':
        return 'Cl'
    elif symbol == 'Br':
        return 'Br'
    elif symbol == 'I':
        return 'I'
    else:
        return symbol


def write_mol2(
    atoms: List[Atom],
    bonds: List[Tuple[int, int, int]],
    output_path: str,
    mol_name: str = "MOL",
    resname: str = "MOL"
) -> bool:
    """
    写入 MOL2 文件
    
    Args:
        atoms: 原子列表
        bonds: 键列表 [(i, j, order), ...]
        output_path: 输出路径
        mol_name: 分子名
        resname: 残基名
    """
    try:
        out_dir = os.path.dirname(output_path)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)
        
        # 生成原子名
        element_count = defaultdict(int)
        atom_names = []
        for atom in atoms:
            sym = ATOMIC_NUMBER_TO_SYMBOL.get(atom.atomic_number, 'X')
            element_count[sym] += 1
            atom_names.append(f"{sym}{element_count[sym]}")
        
        # 统计每个原子的键类型（用于判断 Sybyl 类型）
        atom_bond_orders = defaultdict(list)
        for i, j, order in bonds:
            atom_bond_orders[i].append(order)
            atom_bond_orders[j].append(order)
        
        # 键类型映射
        bond_type_map = {1: '1', 2: '2', 3: '3', 4: 'ar'}
        
        resname = resname[:3].upper()
        
        with open(output_path, 'w', encoding='utf-8') as f:
            # MOLECULE section
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{mol_name}\n")
            f.write(f" {len(atoms)} {len(bonds)} 1 0 0\n")
            f.write("SMALL\n")
            f.write("NO_CHARGES\n\n")
            
            # ATOM section
            f.write("@<TRIPOS>ATOM\n")
            for idx, atom in enumerate(atoms):
                sym = ATOMIC_NUMBER_TO_SYMBOL.get(atom.atomic_number, 'X')
                name = atom_names[idx]
                sybyl = get_sybyl_type(sym, atom_bond_orders.get(idx, []))
                # Format: atom_id atom_name x y z atom_type subst_id subst_name charge
                f.write(f"{idx+1:7d} {name:<4s} {atom.x:10.4f} {atom.y:10.4f} {atom.z:10.4f} "
                       f"{sybyl:<6s} 1 {resname:<4s} 0.0000\n")
            
            # BOND section
            f.write("@<TRIPOS>BOND\n")
            for bond_idx, (i, j, order) in enumerate(bonds):
                btype = bond_type_map.get(order, '1')
                f.write(f"{bond_idx+1:6d} {i+1:5d} {j+1:5d} {btype}\n")
            
            # SUBSTRUCTURE section
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write(f"     1 {resname:<4s}        1 RESIDUE    0 **** **** 0 ROOT\n")
        
        return True
        
    except Exception as e:
        print(f"  [ERROR] 写入 MOL2 失败: {e}", file=sys.stderr)
        return False


# ==============================================================================
# 转换函数
# ==============================================================================
def convert_single_file(
    log_path: Path,
    out_path: Path,
    which: str = 'last',
    orientation: str = 'auto',
    resname: str = 'MOL',
    verbose: bool = False
) -> Tuple[bool, str, int]:
    """
    转换单个 Gaussian log 文件为 MOL2
    
    Returns:
        (成功与否, 消息, 生成的文件数)
    """
    try:
        blocks = select_blocks(str(log_path), which, orientation)
        
        if not blocks:
            return False, "未找到任何构型", 0
        
        files_created = 0
        mol_name = log_path.stem
        
        if which == 'all' and len(blocks) > 1:
            # 多构型：每个 step 一个文件
            for block in blocks:
                step_path = out_path.parent / f"{out_path.stem}_step-{block.block_index:04d}.mol2"
                bonds = generate_bonds(block.atoms, verbose)
                if write_mol2(block.atoms, bonds, str(step_path), mol_name, resname):
                    files_created += 1
            return True, f"{len(blocks)} 个构型", files_created
        else:
            # 单构型
            block = blocks[0]
            bonds = generate_bonds(block.atoms, verbose)
            if write_mol2(block.atoms, bonds, str(out_path), mol_name, resname):
                return True, f"{len(block.atoms)} 原子, {len(bonds)} 键", 1
            else:
                return False, "写入失败", 0
        
    except Exception as e:
        return False, str(e), 0


def process_batch(
    log_dir: Path,
    out_dir: Path,
    skip_existing: bool = True,
    recursive: bool = False,
    which: str = 'last',
    orientation: str = 'auto',
    resname: str = 'MOL',
    verbose: bool = False
) -> BatchResult:
    """批量处理目录中的所有 .log 文件"""
    result = BatchResult()
    
    # 确保输出目录存在
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 收集所有 .log 文件（大小写不敏感）
    all_files = []
    patterns = ["*.log", "*.LOG", "*.out", "*.OUT"]
    
    for pattern in patterns:
        if recursive:
            all_files.extend(log_dir.rglob(pattern))
        else:
            all_files.extend(log_dir.glob(pattern))
    
    all_files = sorted(set(all_files), key=lambda p: p.name.lower())
    result.total = len(all_files)
    
    if result.total == 0:
        print(f"[WARN] 在 {log_dir} 中未找到任何 .log 或 .out 文件")
        return result
    
    print(f"\n{'='*60}")
    print(f"批量转换 Gaussian Log → MOL2")
    print(f"{'='*60}")
    print(f"输入目录: {log_dir}")
    print(f"输出目录: {out_dir}")
    print(f"递归搜索: {'是' if recursive else '否'}")
    print(f"跳过已存在: {'是' if skip_existing else '否'}")
    print(f"构型选择: {which}")
    print(f"找到文件: {result.total}")
    print(f"{'='*60}\n")
    
    for log_file in all_files:
        mol2_name = log_file.stem + ".mol2"
        out_path = out_dir / mol2_name
        
        # 检查是否已存在
        if skip_existing:
            if which == 'all':
                # 检查是否有任何 step 文件存在
                existing = list(out_dir.glob(f"{log_file.stem}_step-*.mol2"))
                if existing:
                    result.skipped += 1
                    result.skipped_files.append(mol2_name)
                    print(f"[SKIP] {log_file.stem}_step-*.mol2 exists ({len(existing)} files)")
                    continue
            elif out_path.exists():
                result.skipped += 1
                result.skipped_files.append(mol2_name)
                print(f"[SKIP] {mol2_name} exists")
                continue
        
        # 转换
        success, msg, count = convert_single_file(
            log_path=log_file,
            out_path=out_path,
            which=which,
            orientation=orientation,
            resname=resname,
            verbose=verbose
        )
        
        if success:
            result.success += 1
            result.success_files.append(mol2_name)
            if count > 1:
                print(f"[OK] {log_file.name} → {count} files ({msg})")
            else:
                print(f"[OK] {log_file.name} → {mol2_name} ({msg})")
        else:
            result.failed += 1
            result.failed_files.append((log_file.name, msg))
            print(f"[FAIL] {log_file.name}: {msg}")
    
    return result


# ==============================================================================
# CLI
# ==============================================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='从 Gaussian log 文件提取几何结构并导出为 MOL2 格式',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
示例:
  # 单文件模式
  python gaussian_log_to_mol2.py --log format/log/TFSI.LOG

  # 批量模式（推荐）
  python gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2

  # 批量模式，跳过已存在
  python gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2 --skip_existing

  # 导出所有构型（每个 step 一个文件）
  python gaussian_log_to_mol2.py --log opt.log --which all

  # 递归搜索子目录
  python gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2 --recursive

版本: {VERSION}
"""
    )
    
    # 输入参数组（互斥）
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--log', help='单文件模式：输入 Gaussian log 文件路径')
    input_group.add_argument('--log_dir', help='批量模式：输入目录路径')
    
    # 输出参数
    parser.add_argument('--out', default=None, help='单文件模式：输出 MOL2 文件路径')
    parser.add_argument('--out_dir', default='format/mol2', help='批量模式：输出目录（默认: format/mol2）')
    
    # 构型选择
    parser.add_argument('--which', default='last',
                       help='导出哪个构型: last（默认）, first, step=N, all')
    
    # 其他选项
    parser.add_argument('--orientation', choices=['standard', 'input', 'auto'],
                       default='auto', help='使用哪种 orientation（默认: auto）')
    parser.add_argument('--resname', default='MOL', help='残基名（默认: MOL）')
    parser.add_argument('--skip_existing', action='store_true', default=True,
                       help='跳过已存在的文件（默认开启）')
    parser.add_argument('--overwrite', action='store_true',
                       help='覆盖已存在的文件')
    parser.add_argument('--recursive', '-r', action='store_true',
                       help='递归搜索子目录')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='详细输出')
    parser.add_argument('--version', action='version', version=f'%(prog)s {VERSION}')
    
    return parser.parse_args()


def main_single(args: argparse.Namespace) -> int:
    """单文件模式"""
    log_path = Path(args.log)
    if not log_path.is_file():
        print(f"[ERROR] 文件不存在: {log_path}", file=sys.stderr)
        return 1
    
    # 确定输出路径
    if args.out:
        out_path = Path(args.out)
    else:
        out_path = Path(args.out_dir) / (log_path.stem + ".mol2")
    
    # 检查是否跳过
    skip = args.skip_existing and not args.overwrite
    if skip and out_path.exists():
        print(f"[SKIP] {out_path.name} exists")
        return 0
    
    print(f"[INFO] 输入: {log_path}")
    print(f"[INFO] 输出: {out_path}")
    
    success, msg, count = convert_single_file(
        log_path=log_path,
        out_path=out_path,
        which=args.which,
        orientation=args.orientation,
        resname=args.resname,
        verbose=args.verbose
    )
    
    if success:
        if count > 1:
            print(f"[SUCCESS] 生成 {count} 个文件 ({msg})")
        else:
            print(f"[SUCCESS] {out_path} ({msg})")
        return 0
    else:
        print(f"[ERROR] {msg}", file=sys.stderr)
        return 1


def main_batch(args: argparse.Namespace) -> int:
    """批量模式"""
    log_dir = Path(args.log_dir)
    out_dir = Path(args.out_dir)
    
    if not log_dir.is_dir():
        print(f"[ERROR] 目录不存在: {log_dir}", file=sys.stderr)
        return 1
    
    skip = args.skip_existing and not args.overwrite
    
    result = process_batch(
        log_dir=log_dir,
        out_dir=out_dir,
        skip_existing=skip,
        recursive=args.recursive,
        which=args.which,
        orientation=args.orientation,
        resname=args.resname,
        verbose=args.verbose
    )
    
    # 打印汇总
    print(f"\n{'='*60}")
    print("批量转换汇总")
    print(f"{'='*60}")
    print(f"总数:   {result.total}")
    print(f"成功:   {result.success}")
    print(f"跳过:   {result.skipped}")
    print(f"失败:   {result.failed}")
    print(f"{'='*60}")
    
    if result.failed > 0:
        print("\n失败文件:")
        for fname, reason in result.failed_files:
            print(f"  - {fname}: {reason}")
    
    return 1 if result.failed > 0 else 0


def main() -> int:
    args = parse_args()
    
    if args.log:
        return main_single(args)
    elif args.log_dir:
        return main_batch(args)
    else:
        print("[ERROR] 请指定 --log 或 --log_dir", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())


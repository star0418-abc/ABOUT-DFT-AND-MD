#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
gaussian_log_to_pdb.py - Gaussian Log 文件转 PDB 格式

从 Gaussian 输出文件 (.log / .out) 中解析几何结构，导出为标准 PDB 格式。
不依赖 OpenBabel 或其他第三方化学库。

功能特点:
  - 支持 Standard orientation 和 Input orientation 解析
  - 流式读取，支持大文件
  - 支持导出单构型或多构型 (multi-model PDB)
  - 支持批量处理目录中的所有 .log 文件
  - 内置元素周期表 (1-118)
  - 严格遵循 PDB 固定列宽格式

用法示例:
  # 单文件模式：导出最后一个构型（默认）
  python scripts/gaussian_log_to_pdb.py --log PTFMA.log

  # 单文件模式：导出第一个构型
  python scripts/gaussian_log_to_pdb.py --log opt.log --which first

  # 单文件模式：导出第 5 个构型
  python scripts/gaussian_log_to_pdb.py --log opt.log --which step=5

  # 单文件模式：导出所有构型为 multi-model PDB
  python scripts/gaussian_log_to_pdb.py --log opt.log --which all --out trajectory.pdb

  # 批量模式：把 format/log 目录下所有 .log 转为 .pdb 到 molecules 目录
  python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules

  # 批量模式：跳过已存在的文件
  python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --skip_existing

  # 批量模式：递归搜索子目录
  python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --recursive

版本: 2.0.0
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import List, Tuple, Optional, Generator, Dict, NamedTuple
from collections import defaultdict
from dataclasses import dataclass



# ==============================================================================
# 内置元素周期表 (Atomic Number -> Element Symbol)
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
    # 扩展到镧系/锕系常见元素
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',  93: 'Np', 94: 'Pu',
    95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg',
    107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Cn',
    113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og',
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
    block_index: int          # 从 1 开始的块索引
    orientation_type: str     # 'standard' 或 'input'
    atoms: List[Atom]         # 原子列表
    line_number: int          # 起始行号（用于调试）


# ==============================================================================
# 解析函数
# ==============================================================================
def iter_orientation_blocks(
    log_path: str,
    orientation_filter: str = 'auto'
) -> Generator[OrientationBlock, None, None]:
    """
    流式解析 Gaussian log 文件，迭代返回 orientation 块。
    
    Args:
        log_path: Gaussian log 文件路径
        orientation_filter: 'auto' | 'standard' | 'input'
            - auto: 优先 Standard orientation，若无则 Input orientation
            - standard: 只返回 Standard orientation
            - input: 只返回 Input orientation
    
    Yields:
        OrientationBlock 对象
    """
    # 正则：匹配 orientation 标记行
    re_standard = re.compile(r'^\s*Standard orientation:', re.IGNORECASE)
    re_input = re.compile(r'^\s*Input orientation:', re.IGNORECASE)
    re_separator = re.compile(r'^\s*-{50,}')  # 分隔线 (至少50个-)
    # 数据行: Center Number, Atomic Number, Atomic Type, X, Y, Z
    # 格式: "    1          6           0       -0.012345    1.234567   -0.987654"
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
            # 检测 orientation 开始
            if re_standard.match(line):
                # 保存之前的块（如果有）
                if in_block and current_atoms:
                    block_index += 1
                    yield OrientationBlock(
                        block_index=block_index,
                        orientation_type=current_type,
                        atoms=current_atoms.copy(),
                        line_number=start_line
                    )
                # 开始新块
                current_type = 'standard'
                current_atoms = []
                in_block = True
                separator_count = 0
                start_line = line_num
                continue
                
            elif re_input.match(line):
                # 只在 filter 允许时处理 Input orientation
                if orientation_filter == 'standard':
                    continue  # 跳过 Input orientation
                    
                # 保存之前的块（如果有）
                if in_block and current_atoms:
                    block_index += 1
                    yield OrientationBlock(
                        block_index=block_index,
                        orientation_type=current_type,
                        atoms=current_atoms.copy(),
                        line_number=start_line
                    )
                # 开始新块
                current_type = 'input'
                current_atoms = []
                in_block = True
                separator_count = 0
                start_line = line_num
                continue
            
            if not in_block:
                continue
                
            # 检测分隔线
            if re_separator.match(line):
                separator_count += 1
                if separator_count >= 3:
                    # 块结束（第三条分隔线）
                    if current_atoms:
                        block_index += 1
                        yield OrientationBlock(
                            block_index=block_index,
                            orientation_type=current_type,
                            atoms=current_atoms.copy(),
                            line_number=start_line
                        )
                    in_block = False
                    current_atoms = []
                continue
            
            # 解析数据行（在第二条分隔线之后）
            if separator_count >= 2:
                m = re_data.match(line)
                if m:
                    atomic_number = int(m.group(2))
                    x = float(m.group(4))
                    y = float(m.group(5))
                    z = float(m.group(6))
                    current_atoms.append(Atom(atomic_number, x, y, z))
    
    # 处理未正常结束的块（log 可能中断）
    if in_block and current_atoms:
        block_index += 1
        yield OrientationBlock(
            block_index=block_index,
            orientation_type=current_type,
            atoms=current_atoms.copy(),
            line_number=start_line
        )


def collect_blocks(
    log_path: str,
    orientation_filter: str = 'auto'
) -> Tuple[List[OrientationBlock], List[OrientationBlock]]:
    """
    收集所有 orientation 块，按类型分类。
    
    Returns:
        (standard_blocks, input_blocks)
    """
    standard_blocks: List[OrientationBlock] = []
    input_blocks: List[OrientationBlock] = []
    
    for block in iter_orientation_blocks(log_path, orientation_filter):
        if block.orientation_type == 'standard':
            standard_blocks.append(block)
        else:
            input_blocks.append(block)
    
    return standard_blocks, input_blocks


def select_blocks(
    log_path: str,
    which: str,
    orientation: str = 'auto'
) -> List[OrientationBlock]:
    """
    根据用户选择返回对应的 orientation 块。
    
    Args:
        log_path: log 文件路径
        which: 'last' | 'first' | 'step=N' | 'all'
        orientation: 'auto' | 'standard' | 'input'
    
    Returns:
        选中的 OrientationBlock 列表
    """
    standard_blocks, input_blocks = collect_blocks(log_path, orientation)
    
    # 根据 orientation 参数决定使用哪组块
    if orientation == 'standard':
        blocks = standard_blocks
        if not blocks:
            raise ValueError(
                f"在 {log_path} 中未找到 Standard orientation 块。\n"
                "提示：尝试使用 --orientation auto 或 --orientation input"
            )
    elif orientation == 'input':
        blocks = input_blocks
        if not blocks:
            raise ValueError(
                f"在 {log_path} 中未找到 Input orientation 块。\n"
                "提示：尝试使用 --orientation auto 或 --orientation standard"
            )
    else:  # auto: 优先 standard，其次 input
        if standard_blocks:
            blocks = standard_blocks
        elif input_blocks:
            blocks = input_blocks
        else:
            blocks = []
    
    if not blocks:
        raise ValueError(
            f"在 {log_path} 中未找到任何 orientation 块。\n"
            "可能的原因：\n"
            "  - 该文件不是 Gaussian 输出文件\n"
            "  - Gaussian 任务异常终止（未完成坐标输出）\n"
            "  - 文件编码问题（尝试用文本编辑器检查）"
        )
    
    # 根据 which 参数选择块
    if which == 'all':
        return blocks
    elif which == 'last':
        return [blocks[-1]]
    elif which == 'first':
        return [blocks[0]]
    elif which.startswith('step='):
        try:
            step_num = int(which.split('=')[1])
        except ValueError:
            raise ValueError(f"无效的 --which 参数: {which}，step=N 中 N 必须是整数")
        
        if step_num < 1:
            raise ValueError(f"step 索引必须 >= 1，收到: {step_num}")
        if step_num > len(blocks):
            raise ValueError(
                f"step={step_num} 超出范围，该文件仅有 {len(blocks)} 个 orientation 块"
            )
        return [blocks[step_num - 1]]
    else:
        raise ValueError(
            f"无效的 --which 参数: {which}\n"
            "有效值: last, first, step=N, all"
        )


# ==============================================================================
# PDB 写出函数
# ==============================================================================
def generate_atom_names(atoms: List[Atom]) -> List[str]:
    """
    为每个原子生成唯一的原子名称。
    格式：元素符号 + 计数，如 C1, C2, H1, H2, ...
    
    PDB 原子名格式要求（4字符）：
    - 1-2字符元素通常左对齐
    - 保证唯一性
    """
    element_counts: Dict[str, int] = defaultdict(int)
    names: List[str] = []
    
    for atom in atoms:
        symbol = ATOMIC_NUMBER_TO_SYMBOL.get(atom.atomic_number, 'X')
        element_counts[symbol] += 1
        count = element_counts[symbol]
        # 生成原子名：最多4字符
        name = f"{symbol}{count}"
        if len(name) > 4:
            # 超长时截断
            name = name[:4]
        names.append(name)
    
    return names


def format_pdb_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str
) -> str:
    """
    格式化单行 PDB HETATM 记录。
    
    PDB 格式规范 (固定列宽):
    COLUMNS        DATA TYPE       FIELD         DEFINITION
    -----------------------------------------------------------------------
     1 -  6        Record name     "HETATM"
     7 - 11        Integer         serial        Atom serial number
    13 - 16        Atom            name          Atom name
    17             Character       altLoc        Alternate location indicator
    18 - 20        Residue name    resName       Residue name
    22             Character       chainID       Chain identifier
    23 - 26        Integer         resSeq        Residue sequence number
    27             AChar           iCode         Insertion code
    31 - 38        Real(8.3)       x             X coordinate (Angstroms)
    39 - 46        Real(8.3)       y             Y coordinate (Angstroms)
    47 - 54        Real(8.3)       z             Z coordinate (Angstroms)
    55 - 60        Real(6.2)       occupancy     Occupancy
    61 - 66        Real(6.2)       tempFactor    Temperature factor
    77 - 78        LString(2)      element       Element symbol (right justified)
    79 - 80        LString(2)      charge        Charge on the atom
    """
    # 原子名对齐规则：
    # - 1-2字符元素：从第14列开始（即在name字段中左对齐）
    # - 但为简单起见，我们将名称左对齐在4字符字段内
    atom_name_formatted = f"{atom_name:<4}"
    
    # 元素符号右对齐（2字符）
    element_formatted = f"{element:>2}"
    
    # 构建完整行
    line = (
        f"HETATM{serial:5d} "                    # 1-11: HETATM + serial
        f"{atom_name_formatted}"                  # 13-16: atom name
        f" "                                      # 17: altLoc
        f"{resname:>3}"                           # 18-20: resName
        f" "                                      # 21: space
        f"{chain:1}"                              # 22: chainID
        f"{resseq:4d}"                            # 23-26: resSeq
        f"    "                                   # 27-30: iCode + spaces
        f"{x:8.3f}{y:8.3f}{z:8.3f}"              # 31-54: coordinates
        f"{1.00:6.2f}"                            # 55-60: occupancy
        f"{0.00:6.2f}"                            # 61-66: tempFactor
        f"          "                             # 67-76: spaces
        f"{element_formatted}"                    # 77-78: element
    )
    
    return line


def write_pdb(
    blocks: List[OrientationBlock],
    output_path: str,
    resname: str = 'UNK',
    chain: str = 'A',
    resseq: int = 1
) -> int:
    """
    将 orientation 块写入 PDB 文件。
    
    Args:
        blocks: OrientationBlock 列表
        output_path: 输出 PDB 文件路径
        resname: 残基名（最多3字符）
        chain: 链 ID（单字符）
        resseq: 残基序号
    
    Returns:
        写入的原子总数
    """
    # 确保输出目录存在
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    
    # 截断参数以符合 PDB 格式
    resname = resname[:3].upper()
    chain = chain[0] if chain else 'A'
    
    total_atoms = 0
    multi_model = len(blocks) > 1
    
    with open(output_path, 'w', encoding='utf-8') as f:
        # 写入头部注释
        f.write(f"REMARK   Generated by gaussian_log_to_pdb.py\n")
        f.write(f"REMARK   Source: {os.path.basename(output_path)}\n")
        f.write(f"REMARK   Blocks: {len(blocks)}\n")
        
        for model_num, block in enumerate(blocks, start=1):
            if multi_model:
                f.write(f"MODEL     {model_num:4d}\n")
            
            atom_names = generate_atom_names(block.atoms)
            
            for serial, (atom, atom_name) in enumerate(zip(block.atoms, atom_names), start=1):
                element = ATOMIC_NUMBER_TO_SYMBOL.get(atom.atomic_number, 'X')
                line = format_pdb_line(
                    serial=serial,
                    atom_name=atom_name,
                    resname=resname,
                    chain=chain,
                    resseq=resseq,
                    x=atom.x,
                    y=atom.y,
                    z=atom.z,
                    element=element
                )
                f.write(line + '\n')
                total_atoms += 1
            
            if multi_model:
                f.write("ENDMDL\n")
        
        f.write("END\n")
    
    return total_atoms


# ==============================================================================
# 批量处理统计
# ==============================================================================
@dataclass
class BatchResult:
    """批量处理结果"""
    total: int = 0
    success: int = 0
    skipped: int = 0
    failed: int = 0
    success_files: List[str] = None
    skipped_files: List[str] = None
    failed_files: List[Tuple[str, str]] = None  # (文件名, 错误原因)
    
    def __post_init__(self):
        if self.success_files is None:
            self.success_files = []
        if self.skipped_files is None:
            self.skipped_files = []
        if self.failed_files is None:
            self.failed_files = []


def convert_single_file(
    log_path: Path,
    out_path: Path,
    which: str = 'last',
    orientation: str = 'auto',
    resname: str = 'UNK',
    chain: str = 'A',
    resseq: int = 1,
    quiet: bool = False
) -> Tuple[bool, str]:
    """
    转换单个 Gaussian log 文件为 PDB
    
    Returns:
        (成功与否, 错误信息或成功消息)
    """
    try:
        # 收集块
        standard_blocks, input_blocks = collect_blocks(str(log_path), orientation)
        
        if not standard_blocks and not input_blocks:
            return False, "未找到任何 orientation 块"
        
        # 选择块
        selected_blocks = select_blocks(str(log_path), which, orientation)
        
        if not selected_blocks:
            return False, "未选中任何构型"
        
        # 写入 PDB
        total_atoms = write_pdb(
            blocks=selected_blocks,
            output_path=str(out_path),
            resname=resname,
            chain=chain,
            resseq=resseq
        )
        
        return True, f"{len(selected_blocks[0].atoms)} 原子"
        
    except Exception as e:
        return False, str(e)


def process_batch(
    log_dir: Path,
    out_dir: Path,
    skip_existing: bool = True,
    recursive: bool = False,
    which: str = 'last',
    orientation: str = 'auto',
    resname: str = 'UNK',
    chain: str = 'A',
    resseq: int = 1,
    quiet: bool = False
) -> BatchResult:
    """
    批量处理目录中的所有 .log 文件
    
    Args:
        log_dir: 输入目录
        out_dir: 输出目录
        skip_existing: 是否跳过已存在的文件
        recursive: 是否递归搜索子目录
        which: 导出哪个构型
        orientation: orientation 类型
        resname: 残基名
        chain: 链 ID
        resseq: 残基序号
        quiet: 安静模式
    
    Returns:
        BatchResult 统计结果
    """
    result = BatchResult()
    
    # 确保输出目录存在
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 收集所有 .log 和 .LOG 文件（大小写不敏感）
    all_files = []
    patterns = ["*.log", "*.LOG", "*.out", "*.OUT"]
    
    for pattern in patterns:
        if recursive:
            all_files.extend(log_dir.rglob(pattern))
        else:
            all_files.extend(log_dir.glob(pattern))
    
    # 去重并按名称排序
    all_files = sorted(set(all_files), key=lambda p: p.name.lower())
    
    result.total = len(all_files)
    
    if result.total == 0:
        if not quiet:
            print(f"[WARN] 在 {log_dir} 中未找到任何 .log 或 .out 文件")
        return result
    
    if not quiet:
        print(f"\n{'='*60}")
        print(f"批量转换 Gaussian Log → PDB")
        print(f"{'='*60}")
        print(f"输入目录: {log_dir}")
        print(f"输出目录: {out_dir}")
        print(f"递归搜索: {'是' if recursive else '否'}")
        print(f"跳过已存在: {'是' if skip_existing else '否'}")
        print(f"找到文件: {result.total}")
        print(f"{'='*60}\n")
    
    for log_file in all_files:
        # 确定输出文件名
        pdb_name = log_file.stem + ".pdb"
        out_path = out_dir / pdb_name
        
        # 检查是否已存在
        if skip_existing and out_path.exists():
            result.skipped += 1
            result.skipped_files.append(pdb_name)
            if not quiet:
                print(f"[SKIP] {pdb_name} exists")
            continue
        
        # 转换
        success, msg = convert_single_file(
            log_path=log_file,
            out_path=out_path,
            which=which,
            orientation=orientation,
            resname=resname,
            chain=chain,
            resseq=resseq,
            quiet=True  # 批量模式下单文件静默
        )
        
        if success:
            result.success += 1
            result.success_files.append(pdb_name)
            if not quiet:
                print(f"[OK] {log_file.name} → {pdb_name} ({msg})")
        else:
            result.failed += 1
            result.failed_files.append((log_file.name, msg))
            if not quiet:
                print(f"[FAIL] {log_file.name}: {msg}")
    
    return result


# ==============================================================================
# CLI 和主函数
# ==============================================================================
def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='从 Gaussian log 文件提取几何结构并导出为 PDB 格式',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
示例:
  # 单文件模式：导出最后一个构型（优化后的最终结构）
  python gaussian_log_to_pdb.py --log opt.log

  # 单文件模式：导出第一个构型
  python gaussian_log_to_pdb.py --log opt.log --which first

  # 单文件模式：导出所有构型为 multi-model PDB
  python gaussian_log_to_pdb.py --log opt.log --which all --out trajectory.pdb

  # 批量模式：把 format/log 目录下所有 .log 转为 .pdb 到 molecules 目录
  python gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules

  # 批量模式：跳过已存在的文件
  python gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --skip_existing

  # 批量模式：递归搜索子目录
  python gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --recursive

版本: {VERSION}
"""
    )
    
    # 输入参数组（互斥）
    input_group = parser.add_mutually_exclusive_group(required=True)
    
    input_group.add_argument(
        '--log',
        help='单文件模式：输入 Gaussian log 文件路径 (.log 或 .out)'
    )
    
    input_group.add_argument(
        '--log_dir',
        default=None,
        help='批量模式：输入目录路径（默认: format/log）'
    )
    
    # 输出参数
    parser.add_argument(
        '--out',
        default=None,
        help='单文件模式：输出 PDB 文件路径（默认：同名 .pdb）'
    )
    
    parser.add_argument(
        '--out_dir',
        default='molecules',
        help='批量模式：输出目录（默认: molecules）'
    )
    
    parser.add_argument(
        '--outdir',
        default=None,
        help='[兼容] 单文件模式的输出目录（与 --out 配合使用）'
    )
    
    # 批量模式选项
    parser.add_argument(
        '--skip_existing',
        action='store_true',
        default=True,
        help='批量模式：跳过已存在的 .pdb 文件（默认: True）'
    )
    
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='批量模式：覆盖已存在的文件（与 --skip_existing 相反）'
    )
    
    parser.add_argument(
        '--recursive', '-r',
        action='store_true',
        help='批量模式：递归搜索子目录中的 .log 文件'
    )
    
    # 通用选项
    parser.add_argument(
        '--which',
        default='last',
        help=(
            '导出哪个构型: '
            'last（默认，最后一个）, first（第一个）, step=N（第N个，从1开始）, all（所有）'
        )
    )
    
    parser.add_argument(
        '--orientation',
        choices=['standard', 'input', 'auto'],
        default='auto',
        help=(
            '使用哪种 orientation: '
            'auto（默认，优先standard，无则input）, standard, input'
        )
    )
    
    parser.add_argument(
        '--resname',
        default='UNK',
        help='PDB 残基名（默认: UNK，最多3字符）'
    )
    
    parser.add_argument(
        '--chain',
        default='A',
        help='PDB 链 ID（默认: A）'
    )
    
    parser.add_argument(
        '--resseq',
        type=int,
        default=1,
        help='PDB 残基序号（默认: 1）'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='安静模式，不输出摘要信息'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {VERSION}'
    )
    
    return parser.parse_args()


def main_single(args: argparse.Namespace) -> int:
    """单文件模式"""
    log_path = os.path.abspath(args.log)
    if not os.path.isfile(log_path):
        print(f"[ERROR] 输入文件不存在: {log_path}", file=sys.stderr)
        return 1
    
    # 确定输出路径
    if args.out:
        out_path = args.out
    else:
        base = os.path.splitext(os.path.basename(log_path))[0]
        out_path = base + '.pdb'
    
    # 处理 --outdir（兼容旧版）
    if args.outdir:
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir, exist_ok=True)
        out_path = os.path.join(args.outdir, os.path.basename(out_path))
    
    out_path = os.path.abspath(out_path)
    
    # 收集统计信息
    try:
        standard_blocks, input_blocks = collect_blocks(log_path, args.orientation)
        
        if not args.quiet:
            print(f"[INFO] 输入文件: {log_path}")
            print(f"[INFO] 找到 {len(standard_blocks)} 个 Standard orientation 块")
            print(f"[INFO] 找到 {len(input_blocks)} 个 Input orientation 块")
    except Exception as e:
        print(f"[ERROR] 解析失败: {e}", file=sys.stderr)
        return 1
    
    # 选择要导出的块
    try:
        selected_blocks = select_blocks(log_path, args.which, args.orientation)
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    
    # 写入 PDB
    try:
        total_atoms = write_pdb(
            blocks=selected_blocks,
            output_path=out_path,
            resname=args.resname,
            chain=args.chain,
            resseq=args.resseq
        )
    except Exception as e:
        print(f"[ERROR] 写入 PDB 失败: {e}", file=sys.stderr)
        return 1
    
    # 输出摘要
    if not args.quiet:
        print(f"[INFO] 导出 {len(selected_blocks)} 个构型")
        print(f"[INFO] 每构型原子数: {len(selected_blocks[0].atoms)}")
        print(f"[INFO] 导出类型: {selected_blocks[0].orientation_type}")
        
        if args.which == 'all':
            print(f"[INFO] Multi-model PDB: 共 {len(selected_blocks)} 个 MODEL")
        elif args.which == 'last':
            print(f"[INFO] 导出了最后一个构型 (块 #{selected_blocks[0].block_index})")
        elif args.which == 'first':
            print(f"[INFO] 导出了第一个构型 (块 #{selected_blocks[0].block_index})")
        elif args.which.startswith('step='):
            print(f"[INFO] 导出了指定构型 (块 #{selected_blocks[0].block_index})")
        
        print(f"[SUCCESS] PDB 已写入: {out_path}")
    
    return 0


def main_batch(args: argparse.Namespace) -> int:
    """批量模式"""
    log_dir = Path(args.log_dir)
    out_dir = Path(args.out_dir)
    
    if not log_dir.is_dir():
        print(f"[ERROR] 输入目录不存在: {log_dir}", file=sys.stderr)
        return 1
    
    # 如果指定了 --overwrite，则不跳过已存在的文件
    skip_existing = args.skip_existing and not args.overwrite
    
    result = process_batch(
        log_dir=log_dir,
        out_dir=out_dir,
        skip_existing=skip_existing,
        recursive=args.recursive,
        which=args.which,
        orientation=args.orientation,
        resname=args.resname,
        chain=args.chain,
        resseq=args.resseq,
        quiet=args.quiet
    )
    
    # 打印汇总
    if not args.quiet:
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
    
    # 有失败则返回非零退出码
    if result.failed > 0:
        return 1
    
    return 0


def main():
    """主函数"""
    args = parse_args()
    
    if args.log:
        # 单文件模式
        sys.exit(main_single(args))
    elif args.log_dir:
        # 批量模式
        sys.exit(main_batch(args))
    else:
        print("[ERROR] 请指定 --log 或 --log_dir", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()


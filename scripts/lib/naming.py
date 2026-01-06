# -*- coding: utf-8 -*-

# Author: Star



"""
naming.py - 统一命名规则模块
============================

提供统一的文件命名和路径生成规则。

命名规则:
- 聚合物输出: P + 原文件名 (EGDA.mol2 → PEGDA.mol2)
- 输出目录: 与输入相同目录
"""

import os
from pathlib import Path
from typing import Optional


def get_polymer_output_path(input_path: str) -> str:
    """
    获取聚合物输出路径
    
    规则: 在原文件名前加 P，保持在同一目录
    
    Args:
        input_path: 输入文件路径 (如 mol2/EGDA.mol2)
    
    Returns:
        输出文件路径 (如 mol2/PEGDA.mol2)
    
    Examples:
        >>> get_polymer_output_path("mol2/EGDA.mol2")
        "mol2/PEGDA.mol2"
        >>> get_polymer_output_path("/path/to/MMA.mol2")
        "/path/to/PMMA.mol2"
    """
    input_path = os.path.abspath(input_path)
    input_dir = os.path.dirname(input_path)
    input_basename = os.path.basename(input_path)
    input_name, input_ext = os.path.splitext(input_basename)
    
    output_name = f"P{input_name}{input_ext}"
    return os.path.join(input_dir, output_name)


def get_polymer_name(monomer_name: str) -> str:
    """
    获取聚合物名称
    
    Args:
        monomer_name: 单体名称 (如 EGDA)
    
    Returns:
        聚合物名称 (如 PEGDA)
    """
    return f"P{monomer_name}"


def get_output_dir(input_path: str) -> str:
    """
    获取输出目录 (与输入相同)
    
    Args:
        input_path: 输入文件路径
    
    Returns:
        输出目录路径
    """
    return os.path.dirname(os.path.abspath(input_path))


def get_timestamped_dir(base_dir: str, prefix: str = "run") -> str:
    """
    获取带时间戳的输出目录
    
    Args:
        base_dir: 基础目录 (如 outputs)
        prefix: 前缀 (如 run)
    
    Returns:
        带时间戳的目录路径 (如 outputs/run_20260101_120000)
    """
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(base_dir, f"{prefix}_{timestamp}")


def ensure_dir(path: str) -> str:
    """
    确保目录存在，不存在则创建
    
    Args:
        path: 目录路径
    
    Returns:
        目录路径
    """
    os.makedirs(path, exist_ok=True)
    return path


def get_project_root() -> Path:
    """获取项目根目录"""
    return Path(__file__).parent.parent.parent


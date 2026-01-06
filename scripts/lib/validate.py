# -*- coding: utf-8 -*-

# Author: Star



"""
validate.py - 统一验证模块
==========================

提供统一的文件存在性、物理量校验等功能。
"""

import os
from typing import List, Optional


class ValidationError(Exception):
    """验证错误"""
    pass


def validate_file_exists(filepath: str, description: str = "文件") -> None:
    """
    验证文件存在
    
    Args:
        filepath: 文件路径
        description: 文件描述 (用于错误信息)
    
    Raises:
        ValidationError: 文件不存在
    """
    if not os.path.exists(filepath):
        raise ValidationError(f"{description}不存在: {filepath}")


def validate_file_nonempty(filepath: str, description: str = "文件") -> int:
    """
    验证文件存在且非空
    
    Args:
        filepath: 文件路径
        description: 文件描述
    
    Returns:
        文件大小 (bytes)
    
    Raises:
        ValidationError: 文件不存在或为空
    """
    validate_file_exists(filepath, description)
    
    size = os.path.getsize(filepath)
    if size == 0:
        raise ValidationError(f"{description}为空: {filepath}")
    
    return size


def validate_output_file(
    filepath: str,
    required_content: Optional[List[str]] = None,
    description: str = "输出文件"
) -> int:
    """
    验证输出文件
    
    Args:
        filepath: 文件路径
        required_content: 必需包含的字符串列表
        description: 文件描述
    
    Returns:
        文件大小 (bytes)
    
    Raises:
        ValidationError: 验证失败
    """
    size = validate_file_nonempty(filepath, description)
    
    if required_content:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
        
        for required in required_content:
            if required not in content:
                raise ValidationError(
                    f"{description}缺少必需内容 '{required}': {filepath}"
                )
    
    return size


def validate_mol2_output(filepath: str) -> int:
    """
    验证 MOL2 输出文件
    
    Returns:
        文件大小 (bytes)
    """
    return validate_output_file(
        filepath,
        required_content=["@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND"],
        description="MOL2 文件"
    )


def validate_pdb_output(filepath: str) -> int:
    """
    验证 PDB 输出文件
    
    Returns:
        文件大小 (bytes)
    """
    return validate_output_file(
        filepath,
        required_content=["ATOM", "END"],
        description="PDB 文件"
    )


def validate_positive(value: float, name: str) -> None:
    """验证正数"""
    if value <= 0:
        raise ValidationError(f"{name} 必须为正数，收到: {value}")


def validate_range(value: float, min_val: float, max_val: float, name: str) -> None:
    """验证范围"""
    if not (min_val <= value <= max_val):
        raise ValidationError(
            f"{name} 必须在 [{min_val}, {max_val}] 范围内，收到: {value}"
        )


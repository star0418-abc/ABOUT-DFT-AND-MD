# -*- coding: utf-8 -*-

# Author: Star



"""
recipe.py - 统一 recipe.yaml 解析模块
======================================

提供统一的 recipe.yaml 读取和解析功能。
所有需要读取配方的脚本都应使用此模块。

用法:
    from lib.recipe import load_recipe, get_oligomer_n
    
    config = load_recipe("config/recipe.yaml")
    n = get_oligomer_n(config, cli_override=5)
"""

import os
from typing import Any, Dict, Optional
from pathlib import Path

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


class RecipeError(Exception):
    """配方解析错误"""
    pass


def load_recipe(recipe_path: str) -> Dict[str, Any]:
    """
    加载 recipe.yaml 配置文件
    
    Args:
        recipe_path: recipe.yaml 文件路径
    
    Returns:
        配置字典
    
    Raises:
        RecipeError: 文件不存在或解析失败
    """
    if not HAS_YAML:
        raise RecipeError("PyYAML 未安装，请运行: pip install pyyaml")
    
    path = Path(recipe_path)
    if not path.exists():
        raise RecipeError(f"配方文件不存在: {recipe_path}")
    
    try:
        with open(path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        if config is None:
            raise RecipeError(f"配方文件为空: {recipe_path}")
        
        return config
    except yaml.YAMLError as e:
        raise RecipeError(f"YAML 解析失败: {e}")


def get_oligomer_n(
    config: Optional[Dict[str, Any]] = None,
    cli_override: Optional[int] = None,
    default: int = 3
) -> int:
    """
    获取聚合度 (oligomer_n)
    
    优先级: CLI > recipe.yaml > 默认值
    
    Args:
        config: recipe 配置字典
        cli_override: CLI 指定的值 (最高优先级)
        default: 默认值 (3)
    
    Returns:
        聚合度
    """
    # 优先级 1: CLI
    if cli_override is not None:
        return cli_override
    
    # 优先级 2: recipe.yaml
    if config:
        # 查找 polymerization.oligomer_n
        polymerization = config.get("polymerization", {})
        if isinstance(polymerization, dict):
            n = polymerization.get("oligomer_n")
            if n is not None:
                return int(n)
        
        # 查找 polymer_matrix[].oligomer_n
        polymer_matrix = config.get("polymer_matrix", [])
        if isinstance(polymer_matrix, list):
            for item in polymer_matrix:
                if isinstance(item, dict) and "oligomer_n" in item:
                    return int(item["oligomer_n"])
    
    # 优先级 3: 默认值
    return default


def get_system_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """获取 system 配置"""
    return config.get("system", {})


def get_packmol_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """获取 packmol 配置"""
    packmol = config.get("packmol", {})
    # 设置默认值
    packmol.setdefault("tolerance_A", 2.0)
    packmol.setdefault("box_scale", 1.5)
    packmol.setdefault("seed", 2025)
    packmol.setdefault("filetype", "pdb")
    packmol.setdefault("output_pdb", "gel.pdb")
    return packmol


def get_htpolynet_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """获取 htpolynet 配置"""
    return config.get("htpolynet", {})


def get_gromacs_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """获取 gromacs 配置"""
    return config.get("gromacs", {})


def get_salt_solution(config: Dict[str, Any]) -> list:
    """获取盐溶液配置"""
    return config.get("salt_solution", [])


def get_polymer_matrix(config: Dict[str, Any]) -> list:
    """获取聚合物基质配置"""
    return config.get("polymer_matrix", [])


def get_monomer_mw(config: Dict[str, Any], monomer_name: str) -> Optional[float]:
    """
    获取单体分子量
    
    Args:
        config: 配置字典
        monomer_name: 单体名称 (如 EGDA, MMA)
    
    Returns:
        分子量 (g/mol) 或 None
    """
    polymerization = config.get("polymerization", {})
    monomer_mw_table = polymerization.get("monomer_mw", {})
    return monomer_mw_table.get(monomer_name)


def get_project_root() -> Path:
    """获取项目根目录"""
    # 假设此文件在 scripts/lib/ 下
    return Path(__file__).parent.parent.parent


def get_default_recipe_path() -> Path:
    """获取默认 recipe.yaml 路径"""
    return get_project_root() / "config" / "recipe.yaml"


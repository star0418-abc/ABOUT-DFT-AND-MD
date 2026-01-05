# -*- coding: utf-8 -*-

# Author: Star



"""
scripts/lib - 公共库模块
========================

提供统一的配置解析、文件读写、验证等功能。

模块:
  - recipe: 统一解析 recipe.yaml
  - io_mol: 分子文件读写与校验
  - naming: 统一命名规则
  - validate: 物理/文件校验
  - logging_utils: 统一日志格式
"""

from .recipe import load_recipe, get_oligomer_n
from .io_mol import read_mol2, write_mol2_strict, validate_mol2_output
from .naming import get_polymer_output_path, get_output_dir
from .validate import validate_file_exists, validate_file_nonempty
from .logging_utils import setup_logger, log_step, log_success, log_error
from .connectivity import check_mol2_connectivity, validate_single_molecule, ConnectivityError

__all__ = [
    'load_recipe', 'get_oligomer_n',
    'read_mol2', 'write_mol2_strict', 'validate_mol2_output',
    'get_polymer_output_path', 'get_output_dir',
    'validate_file_exists', 'validate_file_nonempty',
    'setup_logger', 'log_step', 'log_success', 'log_error',
    'check_mol2_connectivity', 'validate_single_molecule', 'ConnectivityError',
]


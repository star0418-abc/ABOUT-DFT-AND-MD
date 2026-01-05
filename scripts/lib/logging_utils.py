# -*- coding: utf-8 -*-

# Author: Star



"""
logging_utils.py - 统一日志模块
================================

提供统一的日志格式和输出。
"""

import sys
from typing import Optional


# ANSI 颜色代码
class Colors:
    RESET = "\033[0m"
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[0;33m"
    BLUE = "\033[0;34m"
    CYAN = "\033[0;36m"
    BOLD = "\033[1m"


def setup_logger(verbose: bool = False):
    """设置日志级别"""
    # 简单实现，可扩展为完整的 logging 配置
    pass


def log_step(step: str, description: str = "") -> None:
    """
    打印步骤信息
    
    Args:
        step: 步骤编号/名称 (如 "1/4")
        description: 步骤描述
    """
    print(f"\n{Colors.CYAN}[Step {step}]{Colors.RESET} {description}")


def log_success(message: str) -> None:
    """打印成功信息"""
    print(f"{Colors.GREEN}✓{Colors.RESET} {message}")


def log_error(message: str) -> None:
    """打印错误信息"""
    print(f"{Colors.RED}✗{Colors.RESET} {message}", file=sys.stderr)


def log_warning(message: str) -> None:
    """打印警告信息"""
    print(f"{Colors.YELLOW}⚠{Colors.RESET} {message}")


def log_info(message: str) -> None:
    """打印信息"""
    print(f"{Colors.BLUE}ℹ{Colors.RESET} {message}")


def print_header(title: str, width: int = 60) -> None:
    """打印标题头"""
    print("")
    print("=" * width)
    print(f" {title}")
    print("=" * width)


def print_summary(title: str, items: dict, width: int = 60) -> None:
    """
    打印摘要
    
    Args:
        title: 摘要标题
        items: {键: 值} 字典
    """
    print("")
    print("=" * width)
    print(f" {title}")
    print("=" * width)
    for key, value in items.items():
        print(f"  {key}: {value}")
    print("=" * width)


def print_file_info(filepath: str, description: str = "文件") -> None:
    """打印文件信息"""
    import os
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        print(f"  {description}: {filepath} ({size} bytes)")
    else:
        print(f"  {description}: {filepath} (不存在)")


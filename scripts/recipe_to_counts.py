#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Star



"""
recipe_to_counts.py - 配方换算模块
===================================
将用户友好的配方描述（wt%、mol/L）换算为 Packmol 所需的分子数 (count)。

功能：
    - 支持 wt%（质量分数）和 mol/L（摩尔浓度）两种指定方式
    - 自动换算为每个组分的分子数
    - 电荷平衡检查
    - 生成 recipe_resolved.yaml（包含换算后的 count）
    - 生成 summary.json（详细计算过程）

用法：
    python3 scripts/recipe_to_counts.py -c configs/recipe_my.yaml -o outputs/my_run
    python3 scripts/recipe_to_counts.py -c configs/recipe_my.yaml -o outputs/my_run --dry-run

版本: 1.0
"""

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import yaml

# ========== 物理常数 ==========
NA = 6.02214076e23  # Avogadro 常数, mol^-1
CM_TO_A = 1e8       # 1 cm = 1e8 Å
CM3_TO_A3 = 1e24    # 1 cm³ = 1e24 Å³
L_TO_CM3 = 1000     # 1 L = 1000 cm³
NM3_TO_L = 1e-24    # 1 nm³ = 1e-24 L


# ========== 错误类 ==========
class RecipeError(Exception):
    """配方错误"""
    pass


class ValidationError(Exception):
    """校验错误"""
    pass


# ========== 数据结构 ==========
@dataclass
class Component:
    """组分信息"""
    id: str
    name: str
    file: str
    mw_g_mol: float
    charge: int = 0
    
    # 输入（三选一）
    wt_pct: Optional[float] = None
    mol_L: Optional[float] = None
    count: Optional[int] = None
    
    # 计算结果
    count_calc: int = 0
    mass_g: float = 0.0
    mol: float = 0.0
    wt_pct_actual: float = 0.0
    
    # 附加信息
    component_type: str = "component"  # salt_cation, salt_anion, solvent, polymer, crosslinker, initiator
    reactive_atoms: list = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "file": self.file,
            "mw_g_mol": self.mw_g_mol,
            "charge": self.charge,
            "count": self.count_calc,
            "mass_g": f"{self.mass_g:.6e}",
            "mol": f"{self.mol:.6e}",
            "wt_pct_actual": round(self.wt_pct_actual, 4),
            "component_type": self.component_type,
        }


@dataclass
class Salt:
    """盐信息"""
    id: str
    name: str
    cation: Component
    anion: Component
    stoich_cation: int = 1
    stoich_anion: int = 1
    
    # 输入
    mol_L: Optional[float] = None
    wt_pct: Optional[float] = None
    count: Optional[int] = None  # 盐式单元数
    
    # 计算结果
    n_formula: int = 0  # 盐式单元数
    n_cation: int = 0
    n_anion: int = 0
    mass_g: float = 0.0
    
    @property
    def mw_formula(self) -> float:
        """盐式单元分子量"""
        return (self.stoich_cation * self.cation.mw_g_mol + 
                self.stoich_anion * self.anion.mw_g_mol)


# ========== 校验函数 ==========
def validate_file_exists(file_path: str, base_dir: Path) -> None:
    """检查文件是否存在"""
    full_path = base_dir / file_path
    if not full_path.exists():
        raise ValidationError(
            f"文件不存在: {full_path}\n"
            f"  请检查 molecules/ 目录是否包含此文件\n"
            f"  当前目录: {base_dir}"
        )


def validate_mw(mw: float, name: str) -> None:
    """检查分子量是否有效"""
    if mw is None or mw <= 0:
        raise ValidationError(
            f"分子量无效: {name} 的 mw_g_mol = {mw}\n"
            f"  请确保 mw_g_mol 是正数"
        )


def validate_charge_balance(components: list[Component], salts: list[Salt]) -> None:
    """检查电荷平衡"""
    total_charge = 0
    charge_details = []
    
    # 计算盐的电荷
    for salt in salts:
        cation_charge = salt.n_cation * salt.cation.charge
        anion_charge = salt.n_anion * salt.anion.charge
        total_charge += cation_charge + anion_charge
        charge_details.append(f"  {salt.name}: cation({salt.n_cation}×{salt.cation.charge:+d}) + anion({salt.n_anion}×{salt.anion.charge:+d})")
    
    # 计算其他组分的电荷
    for comp in components:
        if comp.charge != 0:
            comp_charge = comp.count_calc * comp.charge
            total_charge += comp_charge
            charge_details.append(f"  {comp.name}: {comp.count_calc}×{comp.charge:+d} = {comp_charge:+d}")
    
    if total_charge != 0:
        raise ValidationError(
            f"电荷不平衡! 净电荷 = {total_charge:+d}\n"
            f"电荷明细:\n" + "\n".join(charge_details) + "\n\n"
            f"修复建议:\n"
            f"  1. 检查盐的化学计量比 (stoichiometry)\n"
            f"  2. 检查离子电荷 (charge) 是否正确\n"
            f"  3. 如果使用带电聚合物，确保有对应的反离子"
        )


def validate_wt_pct_sum(components: list[Component], salts: list[Salt], 
                        target: float = 100.0, tolerance: float = 1.0) -> None:
    """检查质量分数总和"""
    # 只检查使用 wt_pct 指定的组分
    wt_pct_components = []
    total_wt_pct = 0.0
    
    for comp in components:
        if comp.wt_pct is not None:
            wt_pct_components.append(f"  {comp.name}: {comp.wt_pct}%")
            total_wt_pct += comp.wt_pct
    
    for salt in salts:
        if salt.wt_pct is not None:
            wt_pct_components.append(f"  {salt.name}: {salt.wt_pct}%")
            total_wt_pct += salt.wt_pct
    
    if wt_pct_components and abs(total_wt_pct - target) > tolerance:
        raise ValidationError(
            f"质量分数总和不正确!\n"
            f"  指定的 wt% 总和: {total_wt_pct:.2f}%\n"
            f"  目标: {target}% ± {tolerance}%\n\n"
            f"wt% 明细:\n" + "\n".join(wt_pct_components) + "\n\n"
            f"修复建议:\n"
            f"  1. 调整各组分的 wt_pct 使总和接近 100%\n"
            f"  2. 或使用 count 直接指定分子数（跳过 wt% 检查）"
        )


# ========== 换算逻辑 ==========
def calculate_volume_L(system_config: dict, total_mass_g: float, density: float) -> float:
    """计算体系体积 (L)"""
    # 优先使用指定的体积
    if "volume_L" in system_config:
        return system_config["volume_L"]
    
    if "box_size_nm" in system_config:
        box_nm = system_config["box_size_nm"]
        volume_nm3 = box_nm ** 3
        return volume_nm3 * NM3_TO_L
    
    # 否则从质量和密度计算
    volume_cm3 = total_mass_g / density
    return volume_cm3 / L_TO_CM3


def resolve_recipe(config: dict, base_dir: Path) -> tuple[list[Component], list[Salt], dict]:
    """
    解析配方，将 wt%/mol_L 换算为 count
    
    Returns:
        (组分列表, 盐列表, 计算摘要)
    """
    system = config.get("system", {})
    packmol = config.get("packmol", {})
    constraints = config.get("constraints", {})
    
    density = system.get("target_density_g_cm3", 1.0)
    basis_mass = system.get("basis_mass_g", 100.0)
    
    components: list[Component] = []
    salts: list[Salt] = []
    
    # ========== 解析盐溶液 ==========
    salt_solution = config.get("salt_solution", [])
    for item in salt_solution:
        item_type = item.get("type", "")
        
        if item_type == "salt":
            # 解析盐
            cation_data = item.get("cation", {})
            anion_data = item.get("anion", {})
            
            validate_file_exists(cation_data.get("file", ""), base_dir)
            validate_file_exists(anion_data.get("file", ""), base_dir)
            validate_mw(cation_data.get("mw_g_mol"), f"{item['name']} cation")
            validate_mw(anion_data.get("mw_g_mol"), f"{item['name']} anion")
            
            cation = Component(
                id=f"{item['id']}_cation",
                name=cation_data.get("name", "cation"),
                file=cation_data["file"],
                mw_g_mol=cation_data["mw_g_mol"],
                charge=cation_data.get("charge", 1),
                component_type="salt_cation"
            )
            
            anion = Component(
                id=f"{item['id']}_anion",
                name=anion_data.get("name", "anion"),
                file=anion_data["file"],
                mw_g_mol=anion_data["mw_g_mol"],
                charge=anion_data.get("charge", -1),
                component_type="salt_anion"
            )
            
            salt = Salt(
                id=item["id"],
                name=item["name"],
                cation=cation,
                anion=anion,
                stoich_cation=cation_data.get("stoichiometry", 1),
                stoich_anion=anion_data.get("stoichiometry", 1),
                mol_L=item.get("mol_L"),
                wt_pct=item.get("wt_pct"),
                count=item.get("count")
            )
            salts.append(salt)
            
        elif item_type == "solvent" or "file" in item:
            # 解析溶剂或其他组分
            validate_file_exists(item.get("file", ""), base_dir)
            validate_mw(item.get("mw_g_mol"), item.get("name", "unknown"))
            
            comp = Component(
                id=item["id"],
                name=item["name"],
                file=item["file"],
                mw_g_mol=item["mw_g_mol"],
                charge=item.get("charge", 0),
                wt_pct=item.get("wt_pct"),
                mol_L=item.get("mol_L"),
                count=item.get("count"),
                component_type="solvent"
            )
            components.append(comp)
    
    # ========== 解析聚合物基质 ==========
    # 获取全局聚合参数
    polymerization = config.get("polymerization", {})
    global_oligomer_n = polymerization.get("oligomer_n", 3)  # 默认三聚体
    monomer_mw_table = polymerization.get("monomer_mw", {})
    
    polymer_matrix = config.get("polymer_matrix", [])
    for item in polymer_matrix:
        validate_file_exists(item.get("file", ""), base_dir)
        
        # 计算分子量：支持 oligomer_n 自动计算
        # 优先级: 显式 mw_g_mol > oligomer_n × monomer_mw > 报错
        mw_g_mol = item.get("mw_g_mol")
        oligomer_n = item.get("oligomer_n", global_oligomer_n)
        
        if mw_g_mol is None:
            # 尝试从 oligomer_n 计算
            name = item.get("name", "unknown")
            # 尝试从名称推断单体 (PEGDA -> EGDA, PEO -> EO)
            monomer_name = name[1:] if name.startswith("P") and len(name) > 1 else name
            monomer_mw = monomer_mw_table.get(monomer_name)
            
            if monomer_mw:
                mw_g_mol = oligomer_n * monomer_mw
                print(f"  [INFO] {name}: 分子量 = {oligomer_n} × {monomer_mw} = {mw_g_mol:.2f} g/mol")
            else:
                raise RecipeError(
                    f"聚合物 {name} 未指定 mw_g_mol，且无法从 oligomer_n 计算\n"
                    f"  请在 polymer_matrix 中添加 mw_g_mol，或在 polymerization.monomer_mw 中添加 {monomer_name}"
                )
        
        validate_mw(mw_g_mol, item.get("name", "unknown"))
        
        comp = Component(
            id=item["id"],
            name=item["name"],
            file=item["file"],
            mw_g_mol=mw_g_mol,
            charge=item.get("charge", 0),
            wt_pct=item.get("wt_pct"),
            mol_L=item.get("mol_L"),
            count=item.get("count"),
            component_type="polymer",
            reactive_atoms=item.get("reactive_atoms", [])
        )
        components.append(comp)
    
    # ========== 解析交联剂 ==========
    crosslinker = config.get("crosslinker", [])
    for item in crosslinker:
        validate_file_exists(item.get("file", ""), base_dir)
        validate_mw(item.get("mw_g_mol"), item.get("name", "unknown"))
        
        comp = Component(
            id=item["id"],
            name=item["name"],
            file=item["file"],
            mw_g_mol=item["mw_g_mol"],
            charge=item.get("charge", 0),
            wt_pct=item.get("wt_pct"),
            mol_L=item.get("mol_L"),
            count=item.get("count"),
            component_type="crosslinker",
            reactive_atoms=item.get("reactive_atoms", [])
        )
        components.append(comp)
    
    # ========== 解析光引发剂 ==========
    photoinitiator = config.get("photoinitiator", [])
    for item in photoinitiator:
        validate_file_exists(item.get("file", ""), base_dir)
        validate_mw(item.get("mw_g_mol"), item.get("name", "unknown"))
        
        comp = Component(
            id=item["id"],
            name=item["name"],
            file=item["file"],
            mw_g_mol=item["mw_g_mol"],
            charge=item.get("charge", 0),
            wt_pct=item.get("wt_pct"),
            mol_L=item.get("mol_L"),
            count=item.get("count"),
            component_type="initiator"
        )
        components.append(comp)
    
    # ========== 计算分子数 ==========
    # 策略：
    # 1. 如果指定了 count，直接使用
    # 2. 如果指定了 mol_L，用体积计算
    # 3. 如果指定了 wt_pct，用基准质量计算
    # 4. 最后缩放到合理的总分子数（默认 ~1000）
    
    # 计算体积（用于 mol_L 换算）
    volume_L = calculate_volume_L(system, basis_mass, density)
    
    # 目标总分子数（用于缩放）
    target_total = system.get("target_total_molecules", 1000)
    
    # 先处理盐（mol_L 优先）
    for salt in salts:
        if salt.count is not None:
            # 直接使用指定的 count
            salt.n_formula = salt.count
        elif salt.mol_L is not None:
            # 用 mol/L 计算（相对比例）
            mol = salt.mol_L * volume_L
            salt.n_formula = round(mol * NA)
        elif salt.wt_pct is not None:
            # 用 wt% 计算（相对比例）
            mass_g = basis_mass * salt.wt_pct / 100.0
            mol = mass_g / salt.mw_formula
            salt.n_formula = round(mol * NA)
        else:
            raise RecipeError(f"盐 {salt.name} 未指定 count、mol_L 或 wt_pct")
        
        # 确保至少有 1 个（初步）
        salt.n_formula = max(1, salt.n_formula)
        salt.n_cation = salt.n_formula * salt.stoich_cation
        salt.n_anion = salt.n_formula * salt.stoich_anion
        salt.cation.count_calc = salt.n_cation
        salt.anion.count_calc = salt.n_anion
    
    # 处理其他组分
    for comp in components:
        if comp.count is not None:
            # 直接使用指定的 count
            comp.count_calc = comp.count
        elif comp.mol_L is not None:
            # 用 mol/L 计算（相对比例）
            mol = comp.mol_L * volume_L
            comp.count_calc = round(mol * NA)
        elif comp.wt_pct is not None:
            # 用 wt% 计算（相对比例）
            mass_g = basis_mass * comp.wt_pct / 100.0
            mol = mass_g / comp.mw_g_mol
            comp.count_calc = round(mol * NA)
        else:
            raise RecipeError(f"组分 {comp.name} 未指定 count、mol_L 或 wt_pct")
        
        # 确保至少有 1 个
        comp.count_calc = max(1, comp.count_calc)
    
    # ========== 缩放到合理的分子数 ==========
    # 检查是否所有组分都使用了 count（不需要缩放）
    has_direct_count = all(comp.count is not None for comp in components) and \
                       all(salt.count is not None for salt in salts)
    
    if not has_direct_count:
        # 计算当前总分子数
        current_total = sum(comp.count_calc for comp in components) + \
                        sum(salt.n_cation + salt.n_anion for salt in salts)
        
        # 缩放因子
        if current_total > 0:
            scale_factor = target_total / current_total
            
            # 缩放各组分
            for salt in salts:
                salt.n_formula = max(1, round(salt.n_formula * scale_factor))
                salt.n_cation = salt.n_formula * salt.stoich_cation
                salt.n_anion = salt.n_formula * salt.stoich_anion
                salt.cation.count_calc = salt.n_cation
                salt.anion.count_calc = salt.n_anion
            
            for comp in components:
                comp.count_calc = max(1, round(comp.count_calc * scale_factor))
    
    # 重新计算质量和摩尔数
    for salt in salts:
        salt.mass_g = (salt.n_cation * salt.cation.mw_g_mol + 
                       salt.n_anion * salt.anion.mw_g_mol) / NA
        salt.cation.mass_g = salt.n_cation * salt.cation.mw_g_mol / NA
        salt.anion.mass_g = salt.n_anion * salt.anion.mw_g_mol / NA
        salt.cation.mol = salt.n_cation / NA
        salt.anion.mol = salt.n_anion / NA
    
    for comp in components:
        comp.mol = comp.count_calc / NA
        comp.mass_g = comp.mol * comp.mw_g_mol
    
    # ========== 计算总质量和实际质量分数 ==========
    total_mass = sum(comp.mass_g for comp in components) + sum(salt.mass_g for salt in salts)
    
    for comp in components:
        comp.wt_pct_actual = (comp.mass_g / total_mass) * 100.0 if total_mass > 0 else 0.0
    
    for salt in salts:
        salt.cation.wt_pct_actual = (salt.cation.mass_g / total_mass) * 100.0 if total_mass > 0 else 0.0
        salt.anion.wt_pct_actual = (salt.anion.mass_g / total_mass) * 100.0 if total_mass > 0 else 0.0
    
    # ========== 计算盒子尺寸（基于缩放后的实际质量） ==========
    volume_cm3 = total_mass / density
    volume_A3 = volume_cm3 * CM3_TO_A3
    box_L_A = volume_A3 ** (1/3)
    box_L_nm = box_L_A / 10.0
    
    # 重新计算体积（用于摩尔浓度）
    volume_L = volume_cm3 / L_TO_CM3
    
    # 应用 box_scale
    box_scale = packmol.get("box_scale", 1.3)
    packmol_box_L_A = box_L_A * box_scale
    packmol_box_L_nm = box_L_nm * box_scale
    packmol_volume_A3 = packmol_box_L_A ** 3
    packmol_density = total_mass / (packmol_volume_A3 / CM3_TO_A3)
    
    # ========== 校验 ==========
    if constraints.get("enforce_charge_balance", True):
        validate_charge_balance(components, salts)
    
    # 检查 wt% 总和（仅警告，因为可能混用 count 和 wt%）
    try:
        validate_wt_pct_sum(
            components, salts,
            constraints.get("wt_pct_sum_target", 100.0),
            constraints.get("wt_pct_tolerance", 1.0)
        )
    except ValidationError as e:
        print(f"[WARN] {e}")
    
    # ========== 生成摘要 ==========
    summary = {
        "timestamp": datetime.now().isoformat(),
        "system": {
            "name": system.get("name", "unnamed"),
            "target_density_g_cm3": density,
            "target_total_molecules": target_total,
            "actual_total_molecules": sum(comp.count_calc for comp in components) + sum(salt.n_cation + salt.n_anion for salt in salts),
            "volume_L": volume_L,
        },
        "box": {
            "target": {
                "L_A": round(box_L_A, 4),
                "L_nm": round(box_L_nm, 4),
                "volume_A3": round(volume_A3, 4),
            },
            "packmol": {
                "box_scale": box_scale,
                "L_A": round(packmol_box_L_A, 4),
                "L_nm": round(packmol_box_L_nm, 4),
                "volume_A3": round(packmol_volume_A3, 4),
                "initial_density_g_cm3": round(packmol_density, 6),
            },
        },
        "mass": {
            "total_g": f"{total_mass:.6e}",
            "total_Da": f"{total_mass * NA:.6e}",
        },
        "salts": [
            {
                "id": salt.id,
                "name": salt.name,
                "n_formula": salt.n_formula,
                "n_cation": salt.n_cation,
                "n_anion": salt.n_anion,
                "mass_g": f"{salt.mass_g:.6e}",
                "cation": salt.cation.to_dict(),
                "anion": salt.anion.to_dict(),
            }
            for salt in salts
        ],
        "components": [comp.to_dict() for comp in components],
        "packmol": {
            "seed": packmol.get("seed", -1),
            "tolerance_A": packmol.get("tolerance_A", 2.0),
            "output_pdb": packmol.get("output_pdb", "gel.pdb"),
        },
    }
    
    return components, salts, summary


def generate_resolved_yaml(config: dict, components: list[Component], 
                           salts: list[Salt], summary: dict, output_path: Path) -> None:
    """生成解析后的 YAML（包含换算后的 count）
    
    格式与 make_packmol_from_recipe.py 兼容
    """
    # 从盐计算等效浓度
    salt_mol_L = 0.0
    if salts:
        salt = salts[0]
        volume_L = summary['system']['volume_L']
        if volume_L > 0:
            salt_mol_L = salt.n_formula / NA / volume_L
    
    # 生成兼容格式
    resolved = {}
    
    # 注释说明（使用 YAML 注释格式）
    resolved["_comment"] = "由 recipe_to_counts.py 自动生成，所有 count 已根据 wt%/mol_L 换算完成"
    
    # recipe 部分
    resolved["recipe"] = {
        "name": config.get("system", {}).get("name", "resolved"),
        "target_density_g_cm3": config.get("system", {}).get("target_density_g_cm3", 1.0),
        "box_shape": config.get("system", {}).get("box_shape", "cubic"),
        "salt_concentration_mol_L": round(salt_mol_L, 4),  # make_packmol 需要这个字段
        "tolerance_A": config.get("packmol", {}).get("tolerance_A", 2.0),
        "packmol_box_scale": config.get("packmol", {}).get("box_scale", 1.3),
        "packmol_seed": config.get("packmol", {}).get("seed", -1),
        "filetype": "pdb",
        "output_pdb": config.get("packmol", {}).get("output_pdb", "gel.pdb"),
    }
    
    # constraints 部分
    resolved["constraints"] = config.get("constraints", {
        "enforce_mass_balance": True,
        "wt_pct_sum_target": 100.0,
        "wt_pct_tolerance": 1.0,
        "min_count": 1,
        "max_count": 100000,
    })
    
    # salt 部分
    if salts:
        salt = salts[0]  # 目前只支持单盐
        resolved["salt"] = {
            "name": salt.name,
            "cation": {
                "file": salt.cation.file,
                "mw_g_mol": salt.cation.mw_g_mol,
            },
            "anion": {
                "file": salt.anion.file,
                "mw_g_mol": salt.anion.mw_g_mol,
            },
            "stoichiometry": {
                "cation": salt.stoich_cation,
                "anion": salt.stoich_anion,
            },
        }
    
    # components 部分 - 使用 by_count 模式
    resolved["components"] = []
    for comp in components:
        comp_dict = {
            "id": comp.id,
            "name": comp.name,
            "file": comp.file,
            "mw_g_mol": comp.mw_g_mol,
            "mode": "by_count",
            "count": comp.count_calc,
        }
        resolved["components"].append(comp_dict)
    
    # 写入文件（添加头部注释）
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# =============================================================================\n")
        f.write("# 由 recipe_to_counts.py 自动生成\n")
        f.write("# 所有 count 已根据 wt%/mol_L 换算完成\n")
        f.write("# 此文件可直接用于 make_packmol_from_recipe.py\n")
        f.write("# =============================================================================\n\n")
        yaml.dump(resolved, f, allow_unicode=True, default_flow_style=False, sort_keys=False)
    
    print(f"[INFO] 生成 recipe_resolved.yaml: {output_path}")


def generate_summary_json(summary: dict, output_path: Path) -> None:
    """生成摘要 JSON"""
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    
    print(f"[INFO] 生成 summary.json: {output_path}")


# ========== 主函数 ==========
def main():
    parser = argparse.ArgumentParser(
        description="配方换算模块 - 将 wt%/mol_L 换算为分子数 (count)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    # 换算配方并生成输出
    python3 scripts/recipe_to_counts.py -c configs/recipe_my.yaml -o outputs/my_run
    
    # 干运行（仅检验，不生成文件）
    python3 scripts/recipe_to_counts.py -c configs/recipe_my.yaml --dry-run

输出文件:
    outputs/<run>/recipe_resolved.yaml  - 换算后的配方（可直接用于 make_packmol_from_recipe.py）
    outputs/<run>/summary.json          - 详细计算摘要
        """
    )
    
    parser.add_argument(
        "-c", "--config",
        required=True,
        help="配方 YAML 文件路径"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="outputs/resolved",
        help="输出目录 (默认: outputs/resolved)"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="仅校验，不生成文件"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="详细输出"
    )
    
    args = parser.parse_args()
    
    # 读取配置
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"[ERROR] 配置文件不存在: {config_path}")
        return 1
    
    # 推断 base_dir（支持嵌套的 configs 目录）
    base_dir = None
    search_dir = config_path.parent
    for _ in range(5):  # 最多向上搜索 5 级
        if (search_dir / "molecules").exists():
            base_dir = search_dir
            break
        if search_dir.name == "configs":
            candidate = search_dir.parent
            if (candidate / "molecules").exists():
                base_dir = candidate
                break
        if search_dir.parent == search_dir:
            break
        search_dir = search_dir.parent
    
    if base_dir is None:
        base_dir = Path.cwd()
        print(f"[WARN] 无法自动检测项目根目录，使用当前目录: {base_dir}")
    
    print(f"[INFO] 配置文件: {config_path}")
    print(f"[INFO] 基础目录: {base_dir}")
    print("")
    
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    
    try:
        # 解析配方
        components, salts, summary = resolve_recipe(config, base_dir)
        
        # 打印摘要
        print("=" * 60)
        print(" 配方换算结果")
        print("=" * 60)
        print("")
        
        print("盐:")
        for salt in salts:
            salt_wt_pct = (salt.mass_g / float(summary['mass']['total_g'].replace('e', 'E'))) * 100 if float(summary['mass']['total_g'].replace('e', 'E')) > 0 else 0
            print(f"  {salt.name}: {salt.n_formula} 盐式单元 (wt% = {salt_wt_pct:.2f}%)")
            print(f"    阳离子 ({salt.cation.name}): {salt.n_cation}")
            print(f"    阴离子 ({salt.anion.name}): {salt.n_anion}")
        print("")
        
        print("组分:")
        for comp in components:
            print(f"  {comp.name}: {comp.count_calc} 个 (wt% = {comp.wt_pct_actual:.2f}%)")
        print("")
        
        print("盒子:")
        print(f"  目标盒子: {summary['box']['target']['L_nm']:.4f} nm")
        print(f"  Packmol 盒子: {summary['box']['packmol']['L_nm']:.4f} nm (scale={summary['box']['packmol']['box_scale']})")
        print(f"  Packmol 初始密度: {summary['box']['packmol']['initial_density_g_cm3']:.4f} g/cm³")
        print("")
        
        if args.dry_run:
            print("[INFO] 干运行完成，未生成文件")
            return 0
        
        # 创建输出目录
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成文件
        generate_resolved_yaml(config, components, salts, summary, output_dir / "recipe_resolved.yaml")
        generate_summary_json(summary, output_dir / "summary.json")
        
        print("")
        print("=" * 60)
        print(" ✓ 换算完成")
        print("=" * 60)
        print(f"  recipe_resolved.yaml: {output_dir / 'recipe_resolved.yaml'}")
        print(f"  summary.json: {output_dir / 'summary.json'}")
        print("")
        print("下一步:")
        print(f"  python3 scripts/make_packmol_from_recipe.py -c {output_dir / 'recipe_resolved.yaml'} -o {output_dir}")
        
        return 0
        
    except (RecipeError, ValidationError) as e:
        print(f"\n[ERROR] {e}")
        return 1
    except Exception as e:
        print(f"\n[ERROR] 未知错误: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())


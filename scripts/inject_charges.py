# Author: Star


import argparse
import glob
import os
import sys

# ================= 配置区域 =================
MOLECULES_DIR = "molecules"
LOGS_DIR = "gaussian_logs"

# 这里定义哪些分子需要注入 CM5 电荷
# 格式: "Mol2文件名(不带后缀)": "Log文件名关键词"
# 留空则自动尝试在 logs 目录中搜索包含分子名的 log 文件
MOLECULE_MAPPING = {
    "FEC": "FEC",
    "PTFMA": "PTFMA",
    "TFSI": "TFSI",
    "FEP": "FEP",
    # "Li": "Li"  <-- 注意：锂离子不需要注入，请使用手动编写的 ITP
}
# ===========================================

PERIODIC_TABLE = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
    35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe",
}


def parse_mol2_elements(mol2_path):
    """从 MOL2 文件解析元素序列"""
    elements = []
    atom_section = False
    with open(mol2_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped == "@<TRIPOS>ATOM":
                atom_section = True
                continue
            if stripped.startswith("@<TRIPOS>"):
                atom_section = False
                continue
            if atom_section and stripped:
                parts = stripped.split()
                if len(parts) >= 6:
                    atom_type = parts[5]
                    element = atom_type.split(".")[0]
                    element = element.capitalize()
                    elements.append(element)
    return elements


def parse_gaussian_elements(log_file):
    """从 Gaussian log 解析最后一次 Standard/Input orientation 原子序列"""
    try:
        with open(log_file, "r") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"   ❌ 读取 Log 失败: {e}")
        return None

    start_idx = None
    for i in range(len(lines) - 1, -1, -1):
        if "Standard orientation" in lines[i] or "Input orientation" in lines[i]:
            start_idx = i
            break

    if start_idx is None:
        print(f"   ⚠️  警告: 在 {os.path.basename(log_file)} 中未找到坐标块")
        return None

    # 块结构：标题行 + 4 行分隔/表头
    data_start = start_idx + 5
    elements = []
    for line in lines[data_start:]:
        if line.strip().startswith("-----"):
            break
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            atomic_number = int(parts[1])
        except ValueError:
            continue
        element = PERIODIC_TABLE.get(atomic_number)
        if element:
            elements.append(element)
    return elements


def parse_cm5_charges(log_file):
    """从 Gaussian log 解析 CM5 电荷 (读取最后一次出现的电荷)"""
    charges = []
    found_cm5 = False
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            
        # 倒序搜索，确保读取几何优化后的最后一次电荷分布
        start_idx = -1
        for i in range(len(lines) - 1, -1, -1):
            if "CM5 charges" in lines[i] and "Mulliken charges" not in lines[i]:
                start_idx = i + 2 # 跳过标题行
                found_cm5 = True
                break
        
        if not found_cm5:
            print(f"   ⚠️  警告: 在 {os.path.basename(log_file)} 中未找到 'CM5 charges'")
            return None

        # 开始读取
        for i in range(start_idx, len(lines)):
            line = lines[i].strip()
            if "Sum of CM5 charges" in line or line == "":
                break
            
            parts = line.split()
            # Gaussian CM5 格式通常是: Atom No Charge (3列) 或 Atom Charge (2列)
            # 我们取最后一列作为电荷
            if len(parts) >= 2:
                try:
                    charge = float(parts[-1])
                    charges.append(charge)
                except ValueError:
                    continue
                    
        return charges
    except Exception as e:
        print(f"   ❌ 读取 Log 失败: {e}")
        return None

def inject_to_mol2(mol2_path, charges):
    """将电荷写入 Mol2 文件"""
    temp_lines = []
    atom_section = False
    atom_idx = 0
    
    with open(mol2_path, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if line.strip() == "@<TRIPOS>ATOM":
            atom_section = True
            temp_lines.append(line)
            continue
        elif line.strip().startswith("@<TRIPOS>"):
            atom_section = False
            temp_lines.append(line)
            continue
            
        if atom_section:
            parts = line.split()
            if len(parts) >= 9: # 合法的原子行
                if atom_idx >= len(charges):
                    print(f"   ❌ 错误: Mol2 原子数 ({atom_idx+1}+) 超过 Log 电荷数 ({len(charges)})")
                    return False
                
                # 获取新电荷
                new_q = f"{charges[atom_idx]:.6f}"
                
                # 重建行：保留原有内容，只替换最后一列
                # Mol2 格式非常灵活，用空格分隔即可
                # 原始部分: id name x y z type res_id res_name [charge]
                # 我们重新拼接前8列 + 新电荷
                prefix = "\t".join(parts[:-1]) # 使用制表符或空格重新连接前8项
                # 为了保持格式美观，尽量用固定宽度，或者简单地用 split 重组
                
                # 更稳健的方法：
                # line_content = line[:line.rfind(parts[-1])] # 截取到最后一个电荷之前
                # new_line = line_content + new_q + "\n"
                
                # 简单粗暴重组法 (Antechamber 生成的 mol2 通常对齐良好，但 Acpype 读取不挑剔)
                new_line = f"{parts[0]:>7} {parts[1]:<8} {parts[2]:>10} {parts[3]:>10} {parts[4]:>10} {parts[5]:<6} {parts[6]:>4} {parts[7]:<6} {new_q:>10}\n"
                
                temp_lines.append(new_line)
                atom_idx += 1
            else:
                temp_lines.append(line) # 空行或注释
        else:
            temp_lines.append(line) # 非原子部分
            
    if atom_idx != len(charges):
        print(f"   ❌ 错误: 原子数量不匹配! Mol2: {atom_idx}, Log: {len(charges)}")
        print("      请检查 Gaussian 输入文件与 PDB 文件的原子顺序是否一致！")
        return False
        
    # 覆盖原文件
    with open(mol2_path, 'w') as f:
        f.writelines(temp_lines)
    return True

def main():
    parser = argparse.ArgumentParser(description="注入 Gaussian CM5 电荷到 MOL2 文件")
    parser.add_argument("--mol2-dir", default=MOLECULES_DIR, help="MOL2 文件目录")
    parser.add_argument("--log-dir", default=LOGS_DIR, help="Gaussian log 目录")
    parser.add_argument("--strict-mapping", action="store_true", default=True,
                        help="启用原子类型一致性检查 (默认开启)")
    parser.add_argument("--no-strict-mapping", dest="strict_mapping", action="store_false",
                        help="关闭原子类型一致性检查")
    args = parser.parse_args()

    print("🧪 开始注入 CM5 电荷...\n")
    
    for mol_name, log_keyword in MOLECULE_MAPPING.items():
        mol2_path = os.path.join(args.mol2_dir, f"{mol_name}.mol2")
        
        # 1. 检查 Mol2 是否存在
        if not os.path.exists(mol2_path):
            print(f"⏩ 跳过 {mol_name}: 找不到 {mol2_path}")
            continue
            
       # 2. 寻找对应的 Log 文件
        # 尝试 .log (小写)
        log_pattern = os.path.join(args.log_dir, f"*{log_keyword}*.log")
        found_logs = glob.glob(log_pattern)
        
        # 尝试 .LOG (大写) <--- 新增这几行
        if not found_logs:
            log_pattern = os.path.join(args.log_dir, f"*{log_keyword}*.LOG")
            found_logs = glob.glob(log_pattern)

        # 尝试 .out
        if not found_logs:
            log_pattern = os.path.join(args.log_dir, f"*{log_keyword}*.out")
            found_logs = glob.glob(log_pattern)
            
        if not found_logs:
            print(f"⚠️  跳过 {mol_name}: 在 {LOGS_DIR} 未找到包含 '{log_keyword}' 的 log 文件")
            continue
            
        # 默认取找到的第一个 log
        log_file = found_logs[0]
        
        print(f"🔄 处理 {mol_name} ...")
        print(f"   📂 Mol2: {os.path.basename(mol2_path)}")
        print(f"   📄 Log : {os.path.basename(log_file)}")

        if args.strict_mapping:
            mol2_elements = parse_mol2_elements(mol2_path)
            log_elements = parse_gaussian_elements(log_file)
            if log_elements is None:
                print("   ❌ 无法进行原子类型校验，请检查 log 文件")
                continue
            if len(mol2_elements) != len(log_elements):
                print(f"   ❌ 原子数不匹配: Mol2={len(mol2_elements)} Log={len(log_elements)}")
                continue
            if mol2_elements != log_elements:
                print("   ❌ 原子类型顺序不一致，已拒绝注入电荷")
                print(f"      Mol2: {mol2_elements[:10]}{'...' if len(mol2_elements) > 10 else ''}")
                print(f"      Log : {log_elements[:10]}{'...' if len(log_elements) > 10 else ''}")
                continue

        # 3. 解析
        charges = parse_cm5_charges(log_file)
        if not charges:
            continue
            
        # 4. 注入
        if inject_to_mol2(mol2_path, charges):
            print(f"   ✅ 注入成功！")
        else:
            print(f"   ❌ 注入失败。")
        print("-" * 30)

if __name__ == "__main__":
    main()

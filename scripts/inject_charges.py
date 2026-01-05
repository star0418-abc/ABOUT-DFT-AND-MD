# Author: Star


import os
import glob
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
    print("🧪 开始注入 CM5 电荷...\n")
    
    for mol_name, log_keyword in MOLECULE_MAPPING.items():
        mol2_path = os.path.join(MOLECULES_DIR, f"{mol_name}.mol2")
        
        # 1. 检查 Mol2 是否存在
        if not os.path.exists(mol2_path):
            print(f"⏩ 跳过 {mol_name}: 找不到 {mol2_path}")
            continue
            
       # 2. 寻找对应的 Log 文件
        # 尝试 .log (小写)
        log_pattern = os.path.join(LOGS_DIR, f"*{log_keyword}*.log")
        found_logs = glob.glob(log_pattern)
        
        # 尝试 .LOG (大写) <--- 新增这几行
        if not found_logs:
            log_pattern = os.path.join(LOGS_DIR, f"*{log_keyword}*.LOG")
            found_logs = glob.glob(log_pattern)

        # 尝试 .out
        if not found_logs:
            log_pattern = os.path.join(LOGS_DIR, f"*{log_keyword}*.out")
            found_logs = glob.glob(log_pattern)
            
        if not found_logs:
            print(f"⚠️  跳过 {mol_name}: 在 {LOGS_DIR} 未找到包含 '{log_keyword}' 的 log 文件")
            continue
            
        # 默认取找到的第一个 log
        log_file = found_logs[0]
        
        print(f"🔄 处理 {mol_name} ...")
        print(f"   📂 Mol2: {os.path.basename(mol2_path)}")
        print(f"   📄 Log : {os.path.basename(log_file)}")
        
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

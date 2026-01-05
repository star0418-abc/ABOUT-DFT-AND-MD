<<<<<<< HEAD
# ABOUT-DFT-AND-MD
VASP PACKMOL HTPOLYNET GROMACS AND SOMETHING ELSE
  This repository contains code I developed for molecular simulations. Primarily, it uses VASP to calculate the Density of States (DOS) and work function, while enabling an automated workflow from PACKMOL to HTPOLYNET and subsequently to GROMACS. As my research focuses on gel electrolytes, components should be input in the order of salt + polymer matrix + crosslinker + initiator. After inputting the components sequentially, you can obtain your desired results/outputs directly from GROMACS. Overall, the workflow is relatively user-friendly. The code has undergone multiple revisions and updates, with numerous custom scripts integrated to meet my specific research needs.
  Special thanks to Dr. Lu Tian for his excellent MULTIWFN software. Please feel free to raise any questions or issues if encountered.
=======
# 凝胶电解质 Packmol 生成器

配方驱动的分子模拟初始结构生成工具，支持：
- **Packmol** 分子打包
- **HTPolyNet** 聚合物交联网络生成
- **GROMACS** MD 模拟
- **Gaussian log** 文件批量转换为 PDB/MOL2

---

## 目录结构

```
gel_packmol/
├── config/                     # ⭐ 配置文件（唯一权威）
│   └── recipe.yaml            # 统一配置入口
├── configs/                    # 旧配置（MDP 文件仍使用）
│   ├── recipe_my.yaml         # 模板参考
│   ├── recipe_wsgpe_repro.yaml
│   └── mdp/                   # GROMACS MDP 文件
│       ├── em.mdp             # 能量最小化
│       ├── nvt.mdp            # NVT 平衡
│       ├── npt.mdp            # NPT 平衡
│       └── npt_pr.mdp         # NPT 生产
├── format/                     # 格式转换文件
│   ├── log/                   # Gaussian 输出 (.LOG)
│   └── mol2/                  # 转换后的 MOL2 文件
├── molecules/                  # 分子 PDB 文件
│   ├── BF4.pdb, BMIM.pdb      # 离子液体
│   ├── TFSI_TRANS.pdb         # 阴离子
│   ├── PEO.pdb, PMMA.pdb      # 聚合物
│   └── ...
├── htpolynet_out/              # HTPolyNet 输出
│   └── PEGDA/                 # PEGDA 交联网络
├── outputs/                    # 运行输出
├── scripts/                    # ⭐ 所有脚本
│   ├── lib/                   # 公共库模块
│   ├── cli/                   # CLI 入口
│   ├── gaussian_log_to_pdb.py # log → PDB
│   ├── gaussian_log_to_mol2.py# log → MOL2
│   ├── mol2_to_polymer_mol2.py# 单体 → 聚合物
│   ├── make_packmol_from_recipe.py
│   ├── run_workflow.sh        # ⭐ 统一工作流入口
│   ├── run_packmol.sh
│   ├── run_htpolynet.sh
│   └── run_gmx.sh
├── tools/                      # 辅助工具
│   └── build_ptfema.py
├── requirements.txt
└── README.md
```

---

## 🚀 快速开始

### 环境依赖

```bash
# Python 3.8+
pip install pyyaml

# RDKit (聚合物/MOL2 处理必需)
conda install -c conda-forge rdkit

# Packmol
sudo apt install packmol
# 或: conda install -c conda-forge packmol

# HTPolyNet (可选，交联网络)
pip install htpolynet

# GROMACS (可选，MD 模拟)
sudo apt install gromacs
```

### 一键完整流程（推荐）

```bash
cd ~/gel_packmol

# 一键运行: PACKMOL → HTPOLYNET → GROMACS
bash scripts/run_workflow.sh

# 显式指定配置文件
bash scripts/run_workflow.sh -c config/recipe.yaml

# 跳过某些阶段
bash scripts/run_workflow.sh --skip-htpolynet --skip-gromacs
```

---

## 📋 常用命令速查

### 1. Gaussian log → PDB

```bash
# ⭐ 批量转换（推荐）
python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules

# 批量转换，跳过已存在
python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --skip_existing

# 批量转换，递归搜索子目录
python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --recursive

# 单文件模式
python scripts/gaussian_log_to_pdb.py --log format/log/TFSI_TRANS.LOG --out molecules/TFSI.pdb

# 导出所有构型（multi-model PDB）
python scripts/gaussian_log_to_pdb.py --log opt.log --which all --out trajectory.pdb
```

**输入**: `format/log/*.LOG` (或 `.log`, `.out`)  
**输出**: `molecules/*.pdb`  
**关键参数**: `--skip_existing`, `--recursive`, `--which [last|first|step=N|all]`

### 2. Gaussian log → MOL2

```bash
# ⭐ 批量转换（推荐）
python scripts/gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2

# 批量转换，跳过已存在
python scripts/gaussian_log_to_mol2.py --log_dir format/log --out_dir format/mol2 --skip_existing

# 单文件模式
python scripts/gaussian_log_to_mol2.py --log format/log/TFSI_TRANS.LOG --out format/mol2/TFSI.mol2

# 导出所有构型（每个 step 一个文件）
python scripts/gaussian_log_to_mol2.py --log opt.log --which all
```

**输入**: `format/log/*.LOG`  
**输出**: `format/mol2/*.mol2`  
**关键参数**: `--skip_existing`, `--which [last|first|step=N|all]`, `--resname`

### 3. 单体 MOL2 → 聚合物 MOL2

```bash
# ⭐ 批量处理四个单体（推荐）
python scripts/process_monomers_batch.py

# 指定单体和聚合度
python scripts/process_monomers_batch.py --names EGDA MMA TFEMA VA --n_repeat 3

# 单文件模式（默认 n=3）
python scripts/mol2_to_polymer_mol2.py --input mol2/MMA.mol2

# 使用 recipe.yaml 配置
python scripts/mol2_to_polymer_mol2.py --input mol2/EGDA.mol2 --recipe config/recipe.yaml

# CLI 覆盖聚合度
python scripts/mol2_to_polymer_mol2.py --input mol2/MMA.mol2 --n_repeat 5

# EO 开环聚合
python scripts/mol2_to_polymer_mol2.py --input mol2/EO.mol2 --n_repeat 3 -v
```

**输入**: `mol2/X.mol2` (或任意路径)  
**输出**: 同目录 `PX.mol2`  
**支持的单体**: EO→PEO, VA→PVA, TFEMA→PTFEMA, MMA→PMMA, AM→PAM, EGDA→PEGDA

### 4. 从 recipe.yaml 生成 Packmol 输入

```bash
# 生成 Packmol 输入 + 运行 Packmol
python scripts/make_packmol_from_recipe.py -c config/recipe.yaml -o outputs

# 仅计算，不生成文件
python scripts/make_packmol_from_recipe.py -c config/recipe.yaml --dry-run

# 运行 Packmol（使用脚本）
bash scripts/run_packmol.sh -c config/recipe.yaml -o outputs/my_run
```

**输入**: `config/recipe.yaml`  
**输出**: `outputs/gel.inp`, `outputs/summary.json`, `outputs/gel.pdb`

### 5. HTPolyNet 交联网络

```bash
# 生成 PEGDA 交联网络
python scripts/build_pegda_network.py

# 自定义参数
python scripts/build_pegda_network.py --n_monomers 300 --conversion 0.9

# 运行 HTPolyNet 官方示例（验证环境）
bash scripts/run_htpolynet.sh --example 0 -o outputs/htpolynet_test
```

**输出**: `htpolynet_out/PEGDA/` (pdb/gro/top)

### 6. GROMACS 模拟

```bash
# 运行 GROMACS (EM → NVT → NPT)
bash scripts/run_gmx.sh -i outputs/htpolynet/PEGDA -o outputs/gmx
```

---

## 📁 scripts/ 目录详解

| 脚本 | 用途 | 输入 | 输出 |
|------|------|------|------|
| `gaussian_log_to_pdb.py` | Gaussian log 批量转 PDB | `format/log/*.LOG` | `molecules/*.pdb` |
| `gaussian_log_to_mol2.py` | Gaussian log 批量转 MOL2 | `format/log/*.LOG` | `format/mol2/*.mol2` |
| `mol2_to_polymer_mol2.py` | 单体 MOL2 → 聚合物 MOL2 | `mol2/X.mol2` | `mol2/PX.mol2` |
| `process_monomers_batch.py` | 批量处理四个单体 | 配置/CLI | `mol2/P*.mol2` |
| `make_packmol_from_recipe.py` | 生成 Packmol 输入 | `config/recipe.yaml` | `outputs/gel.inp` |
| `recipe_to_counts.py` | 配方换算 (wt%/mol/L → count) | `config/recipe.yaml` | `recipe_resolved.yaml` |
| `build_pegda_network.py` | PEGDA 交联网络 (HTPolyNet) | - | `htpolynet_out/PEGDA/` |
| `build_pva_trimer.py` | 生成 PVA 三聚体 | `mol2/VA.mol2` | `mol2/PVA.mol2` |
| `generate_topology.py` | 生成 GROMACS 拓扑 (GAFF) | PDB + recipe | `topol.top` |
| `generate_topology_gaff.py` | 生成 GROMACS 拓扑 (GAFF2) | PDB + recipe | `topol.top` |
| `inject_charges.py` | 注入 CM5 电荷到 mol2 | log + mol2 | mol2 (覆盖) |
| `selfcheck.py` | 整体自检 | - | 测试报告 |
| `selfcheck_oligomer.py` | 寡聚体生成自检 | - | 测试报告 |
| `sanitize_pdb_for_ligpargen.py` | 清理 PDB for LigParGen | PDB | PDB |

### Shell 脚本

| 脚本 | 用途 | 推荐用法 |
|------|------|----------|
| `run_workflow.sh` | ⭐ 统一工作流入口 | `bash scripts/run_workflow.sh` |
| `run_packmol.sh` | Packmol 驱动 | `bash scripts/run_packmol.sh -c config/recipe.yaml -o outputs` |
| `run_htpolynet.sh` | HTPolyNet 驱动 | `bash scripts/run_htpolynet.sh --example 0 -o outputs/htpolynet` |
| `run_gmx.sh` | GROMACS 驱动 | `bash scripts/run_gmx.sh -i outputs/htpolynet -o outputs/gmx` |
| `run_full.sh` | 完整流程 (旧) | 已弃用，请用 `run_workflow.sh` |
| `run_all_smoke.sh` | Smoke 测试 (旧) | 已弃用 |
| `open_in_windows.sh` | 在 Windows 打开目录 | `bash scripts/open_in_windows.sh outputs/` |

### lib/ 公共库

| 模块 | 用途 |
|------|------|
| `recipe.py` | 统一 recipe.yaml 解析 |
| `io_mol.py` | 分子文件读写与校验 |
| `naming.py` | 统一命名规则 |
| `validate.py` | 文件/物理量校验 |
| `logging_utils.py` | 统一日志格式 |
| `connectivity.py` | 分子连通性检查 |

---

## ⚠️ 常见问题 (Troubleshooting)

### 1. 文件找不到 (FileNotFoundError)

```
FileNotFoundError: [Errno 2] No such file or directory: 'format/log/XXX.log'
```

**原因**: 脚本默认路径与实际文件位置不符  
**解决**: 检查文件是否存在，使用 `--log_dir` 或 `--log` 显式指定路径

```bash
# 检查文件
ls format/log/*.LOG

# 显式指定
python scripts/gaussian_log_to_pdb.py --log_dir /path/to/your/logs
```

### 2. 输出目录不存在

```
FileNotFoundError: [Errno 2] No such file or directory: 'molecules/XXX.pdb'
```

**原因**: 输出目录不存在  
**解决**: 脚本会自动创建目录，但如果权限不足则失败

```bash
# 手动创建
mkdir -p molecules format/mol2 outputs
```

### 3. 同名文件覆盖风险

**问题**: 批量转换时可能覆盖已有文件  
**解决**: 使用 `--skip_existing`（大多数脚本默认开启）

```bash
# 安全模式（推荐）
python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --skip_existing

# 强制覆盖
python scripts/gaussian_log_to_pdb.py --log_dir format/log --out_dir molecules --overwrite
```

### 4. 解析不到 orientation

```
[FAIL] XXX.LOG: 未找到任何 orientation 块
```

**原因**: Gaussian log 文件格式异常或不包含几何优化  
**解决**: 检查 log 文件内容，确认包含 "Standard orientation" 或 "Input orientation"

```bash
grep -c "orientation" format/log/XXX.LOG
```

### 5. 空目录运行

```
============================================================
批量转换汇总
============================================================
总数:   0
成功:   0
跳过:   0
失败:   0
============================================================
```

**解释**: 目录中没有匹配的 .log/.LOG 文件，脚本正常退出  
**行为**: 退出码为 0（成功），只是没有文件需要处理

### 6. RDKit 未安装

```
错误: RDKit 未安装
```

**解决**: 安装 RDKit

```bash
conda install -c conda-forge rdkit
```

### 7. Packmol Density Trap

```
ERROR: Packmol could not find a solution!
```

**原因**: 盒子太小，分子无法放置  
**解决**: 在 recipe.yaml 中增大 `packmol_box_scale`

```yaml
recipe:
  packmol_box_scale: 1.5  # 从 1.3 增到 1.5
```

### 8. HTPolyNet 无法识别反应位点

```
ERROR: No reaction sites found
```

**原因**: PDB 原子命名与 HTPolyNet cfg 不匹配  
**解决**: 检查 PDB 原子名，确保双键碳命名正确

---

## 配置文件说明

### config/recipe.yaml（唯一权威配置）

```yaml
# 体系参数
recipe:
  name: "LiTFSI-PEGDA-gel"
  temperature_K: 298
  target_density_g_cm3: 1.15
  salt_concentration_mol_L: 1.0
  packmol_box_scale: 1.3       # 盒子膨胀系数

# 组分定义
components:
  - id: "polymer"
    name: "PEGDA"
    file: "molecules/PEGDA.pdb"
    count: 10
    mw_g_mol: 700.0

  - id: "solvent"
    name: "EC"
    file: "molecules/EC.pdb"
    mode: "by_wt_pct"
    wt_pct: 30.0

# 盐定义
salt:
  name: "LiTFSI"
  stoichiometry: {cation: 1, anion: 1}
  cation: {file: "molecules/Li.pdb", mw_g_mol: 6.94}
  anion: {file: "molecules/TFSI.pdb", mw_g_mol: 280.13}
```

---

## 自检验证

```bash
cd ~/gel_packmol

# Python 语法检查
python -m py_compile scripts/gaussian_log_to_pdb.py
python -m py_compile scripts/gaussian_log_to_mol2.py
python -m py_compile scripts/mol2_to_polymer_mol2.py

# 功能自检
python scripts/selfcheck.py
python scripts/selfcheck_oligomer.py

# 寡聚体生成自检
python scripts/mol2_to_polymer_mol2.py --selfcheck
```

---

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request!
>>>>>>> e5c09c9 (init)

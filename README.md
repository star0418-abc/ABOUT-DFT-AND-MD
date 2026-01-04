# ABOUT-DFT-AND-MD
VASP PACKMOL HTPOLYNET GROMACS AND SOMETHING ELSE
  This repository contains code I developed for molecular simulations. Primarily, it uses VASP to calculate the Density of States (DOS) and work function, while enabling an automated workflow from PACKMOL to HTPOLYNET and subsequently to GROMACS. As my research focuses on gel electrolytes, components should be input in the order of salt + polymer matrix + crosslinker + initiator. After inputting the components sequentially, you can obtain your desired results/outputs directly from GROMACS. Overall, the workflow is relatively user-friendly. The code has undergone multiple revisions and updates, with numerous custom scripts integrated to meet my specific research needs.
  Special thanks to Dr. Lu Tian for his excellent MULTIWFN software. Please feel free to raise any questions or issues if encountered.
# VASP Scripts - 凝胶电解质计算工作流

针对 VASP (oneAPI + Intel MPI) 环境的计算辅助脚本集合，适用于凝胶电解质、AIMD 等计算任务。

## 文件列表

| 文件 | 用途 |
|------|------|
| `vasp_env.sh` | VASP 运行环境配置（oneAPI、自检、线程控制） |
| `run_vasp.sh` | VASP 运行脚本（备份、日志、续算、核数检查） |
| `check_vasp.sh` | 检查计算状态（完成标志、能量、费米能级） |
| `aimd_watch.sh` | AIMD 实时监控（温度、离子步、能量） |
| `aimd_msd.py` | MSD 计算与扩散系数拟合（需要 numpy） |
| `aimd_post.py` | AIMD 热力学数据后处理（E0、T、F 导出 CSV） |
| `clean_vasp.sh` | 安全清理大文件（WAVECAR、CHGCAR 等） |
| `recipe.yaml` | 配方定义文件示例（8 类组分 + 模拟条件） |
| `recipe_validate.py` | 配方验证工具（校验 wt% 总和、温度、格式） |
| `recipe_to_counts.py` | 配方换算工具（wt% → 分子/原子数） |
| `make_incar_aimd.py` | AIMD INCAR 生成器（核心） |
| `aimd_setup.sh` | AIMD 一键设置脚本 |
| `setup_electronic.py` | 电子性质输入生成（功函数/DOS） |
| `analyze_electronic.py` | 电子性质后处理（功函数/DOS） |
| `setup_aimd_ase.py` | 从大体系切割 AIMD cluster |

## 快速开始

### 1. 安装

```bash
cd ~/vasp_scripts
chmod +x *.sh *.py
source ~/.bashrc
```

### 2. 依赖安装

```bash
pip install numpy pyyaml ase matplotlib
```

---

## 环境配置 (vasp_env.sh)

### 功能

- 加载 oneAPI 环境（失败时明确报错，不静默）
- 自检 VASP 可执行文件（vasp_std/gam/ncl）
- 自检 mpirun 可用性
- 检查 vasp_std 动态库依赖
- 设置默认纯 MPI 模式（OMP_NUM_THREADS=1）

### 用法

```bash
source ~/vasp_scripts/vasp_env.sh
```

### 环境变量

| 变量 | 说明 |
|------|------|
| `VASP_BIN` | VASP 可执行文件目录 |
| `VASP_PP_PATH` | POTCAR 路径（用于自动拼接脚本） |
| `OMP_NUM_THREADS` | OpenMP 线程数（默认 1） |
| `MKL_NUM_THREADS` | MKL 线程数（默认 1） |
| `I_MPI_ADJUST_REDUCE` | Intel MPI 参数（不强制覆盖，尊重用户设置） |

---

## VASP 运行脚本 (run_vasp.sh)

### 功能

- ✅ WSL 核数检查（自动下调或报错）
- ✅ 磁盘空间预检查
- ✅ AIMD 续算支持（CONTCAR → POSCAR）
- ✅ stdout/stderr 分离输出
- ✅ 运行耗时统计
- ✅ 自动备份输入文件
- ✅ 归档旧输出文件
- ✅ 并行参数提示（NCORE/KPAR）

### 用法

```bash
# 基本运行
NP=16 EXE=vasp_std run_vasp.sh

# 续算模式
RESUME=1 NP=16 run_vasp.sh

# 严格核数检查（超限报错）
STRICT_NP=1 NP=32 run_vasp.sh

# 强制忽略磁盘检查
FORCE_DISK=1 run_vasp.sh
```

### 环境变量

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `NP` | 8 | MPI 进程数 |
| `EXE` | vasp_std | 可执行文件 |
| `OUT` | vasp.out | stdout 文件 |
| `ERR` | vasp.err | stderr 文件 |
| `RESUME` | 0 | 续算模式：1=自动 cp CONTCAR→POSCAR |
| `STRICT_NP` | 0 | 严格核数：1=超限报错，0=自动下调 |
| `RESERVE_CORES` | 2 | WSL 预留核数 |
| `MIN_FREE_GB` | 20 | 最小磁盘空间 (GB) |
| `FORCE_DISK` | 1 | 忽略磁盘检查 |

### 输出文件

- `vasp.out` - 标准输出
- `vasp.err` - 标准错误
- `run.log` - 运行日志（含耗时、状态）
- `snapshots/<timestamp>/` - 输入文件备份
- `old/` - 旧输出文件归档

---

## 配方层 (Recipe Layer)

### 配方组分（8 类固定顺序）

| 序号 | 类别 | 中文名 | 必需性 |
|------|------|--------|--------|
| 1 | `salt_solution` | 盐溶液 | 主要 |
| 2 | `polymer_matrix` | 聚合物基质 | 主要 |
| 3 | `crosslinker` | 交联剂 | 主要 |
| 4 | `photoinitiator` | 引发剂 | 主要 |
| 5 | `plasticizer_solvent` | 增塑剂/溶剂 | 可选 |
| 6 | `functional_monomer` | 功能单体 | 可选 |
| 7 | `stabilizer` | 稳定剂 | 可选 |
| 8 | `functional_filler` | 功能填料 | 可选 |

### recipe.yaml 结构

```yaml
# 模拟条件
simulation:
  mode: aimd
  temperature_C: 60       # 摄氏度 → 自动转 K
  dt_fs: 1.0              # 时间步长 = POTIM
  nsteps: 10000           # 总步数 = NSW
  ensemble: nvt
  thermostat: langevin
  gamma_1ps: 10.0         # Langevin 摩擦系数

  # AIMD 稳定性参数
  isym: 0                 # AIMD 必须关闭对称性
  maxmix: 40              # 电荷密度混合历史，建议 40-80

  # 建模参数（非 VASP 参数）
  density_g_cm3: 1.25     # 体系密度，用于建盒子
  builder: none           # 建模工具: none/packmol
  target_atoms: 200       # AIMD 目标原子数
  allow_drop_low_fraction_components: true

# 组分（wt%）
salt_solution:
  - name: "LiTFSI（双三氟甲磺酰亚胺锂）"
    wt_pct: 15.0
    kind: salt
    mw_g_mol: 287.09
    atoms_per_entity: 15

polymer_matrix:
  - name: "PEGDA（聚乙二醇二丙烯酸酯）"
    wt_pct: 40.0
    kind: polymer
    ...

# 可选组分为空
plasticizer_solvent: []
functional_monomer: []
stabilizer: []
functional_filler: []
```

### 命名规范

> ⚠️ **重要**：`name` 字段必须包含 **缩写（中文全称）**
> 
> 例如：`"LiTFSI（双三氟甲磺酰亚胺锂）"`、`"PEGDA（聚乙二醇二丙烯酸酯）"`

### 关于 density 和 target_atoms

- `density_g_cm3` 仅用于建盒子/Packmol/recipe_to_counts，非 VASP 参数
- `target_atoms` 用于验证/提示，AIMD 代表性小胞可能无法严格保持实验 wt%
- `allow_drop_low_fraction_components: true` 允许小体系中 optional 组分 count=0

---

## 从大体系切割 AIMD Cluster（ASE）

### 概述

`setup_aimd_ase.py` 用于从大体系结构（GROMACS/Packmol 输出）中切割出一个可用于 VASP AIMD 的局部量子区域。

**典型场景**：凝胶电解质体系有数千原子，AIMD 只能算几百原子，需要切割一个以目标离子为中心的小 cluster。

### 基本用法

```bash
# 以 Li 原子为中心，半径 8 Å，温度 350 K
python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom Li --radius 8 --temp 350 --outdir aimd_Li8A

# 使用原子索引（0-based）
python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom 100 --radius 10 --outdir aimd_idx100

# 保留完整分子（避免切断分子）
python3 setup_aimd_ase.py --src system.pdb --center_atom Li --radius 8 --selection molecule --outdir aimd_mol
```

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--src` | 必填 | 输入结构文件（.pdb/.gro/.xyz/.cif） |
| `--center_atom` | 必填 | 中心原子（索引或元素符号如 Li） |
| `--radius` | 8.0 | 切割半径 Å |
| `--selection` | sphere | 选择模式: sphere / molecule |
| `--vacuum` | 20 | 真空层厚度 Å |
| `--temp` | 350 | AIMD 温度 K |
| `--steps` | 2000 | AIMD 步数 |
| `--potim` | 2.0 | 时间步长 fs（含 H 建议 0.5-1.0） |
| `--kpoints` | "1 1 1" | K 点网格（Gamma-only） |
| `--ncore` | - | NCORE 并行参数 |
| `--outdir` | aimd_cluster | 输出目录 |
| `--overwrite` | False | 覆盖已存在目录 |
| `--one_based` | False | 原子索引按 1-based 解释 |

### 选择模式

- **sphere**：纯半径选择，最快
- **molecule**：半径内的原子所属分子整体保留，避免切断分子（需要结构文件包含分子信息，如 .pdb 的 residue）

> **建议**：凝胶电解质体系使用 `--selection molecule` 更物理合理

### 输出文件

```
aimd_Li8A/
├── POSCAR              # VASP 结构文件
├── INCAR               # AIMD 参数
├── KPOINTS             # K 点（Gamma-only）
├── POTCAR              # 赝势（如 VASP_PP_PATH 已设置）
├── cluster_visual.xyz  # 可视化文件（OVITO/VMD）
└── selected_indices.txt # 选中的原子索引
```

### 电荷检查

脚本会自动估算 cluster 总电荷（基于离子电荷表）：
- 若 |电荷| ≥ 1，会打印警告
- 可能需要调整 VASP NELECT 或重新选择半径/中心

### 完整工作流

```bash
# 1. 从 GROMACS 输出切割 cluster
python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom Li --radius 8 \
    --selection molecule --vacuum 20 --temp 350 --steps 2000 --outdir aimd_Li8A

# 2. 检查 cluster（可选）
# 用 OVITO 打开 aimd_Li8A/cluster_visual.xyz

# 3. 准备 POTCAR（如未自动生成）
export VASP_PP_PATH=/path/to/potentials
# 或手动拼接

# 4. 运行 AIMD
cd aimd_Li8A
NP=16 EXE=vasp_std run_vasp.sh

# 5. 监控
aimd_watch.sh

# 6. 后处理
python3 aimd_post.py
python3 aimd_msd.py --specie Li --dt_fs 2.0
```

### 注意事项

- **.gro 文件问题**：ASE 对 .gro 支持有限，建议先转换：
  ```bash
  gmx trjconv -f input.gro -o output.pdb
  ```

- **含 H 原子**：时间步长建议 0.5-1.0 fs（默认 2.0 可能过大）

- **PBC 处理**：脚本自动使用最小镜像距离（MIC）选择原子

- **molecule 模式**：需要结构文件包含分子信息（如 .pdb 的 residue），否则回退到 sphere

---

## AIMD 工作流

### 完整流程

```bash
# 1. 准备配方
cp ~/vasp_scripts/recipe.yaml ./
vim recipe.yaml  # 修改温度、组分

# 2. 验证配方
python3 recipe_validate.py

# 3. 生成 INCAR
python3 make_incar_aimd.py --out INCAR
# 或使用一键设置
aimd_setup.sh

# 4. 运行 AIMD
NP=16 EXE=vasp_std run_vasp.sh

# 5. 监控
aimd_watch.sh
# 或
tail -f vasp.out vasp.err

# 6. 续算（如需要）
RESUME=1 NP=16 run_vasp.sh

# 7. 后处理
python3 aimd_post.py
python3 aimd_msd.py --specie Li --dt_fs 1.0
```

### make_incar_aimd.py 功能

- 从 recipe.yaml 读取模拟条件
- 温度自动转换：°C → K
- 强制设置 `ISYM = 0`（AIMD 必须）
- 设置 `MAXMIX`（默认 40）
- 添加 `LASPH = .TRUE.` 和 `ADDGRID = .TRUE.`
- 支持 INCAR.base 继承

### 生成的 INCAR 关键参数

| 参数 | 值 | 说明 |
|------|-----|------|
| IBRION | 0 | MD 模式 |
| NSW | nsteps | 离子步数 |
| POTIM | dt_fs | 时间步长 |
| TEBEG/TEEND | T_K | 温度 (K) |
| ISYM | 0 | 关闭对称性（强制） |
| MAXMIX | 40 | 电荷混合历史 |
| MDALGO | 2/3 | Langevin/Nosé-Hoover |
| LWAVE | .FALSE. | 不写 WAVECAR |
| LCHARG | .FALSE. | 不写 CHGCAR |

---

## 电子性质计算（功函数/DOS）

### 概述

支持两种电子性质计算：
- **功函数 (Work Function)**: Φ = V_vac − E_F，从 LOCPOT + OUTCAR 计算
- **DOS/PDOS**: 态密度，使用两步法（SCF → CHGCAR → NSCF）

### 功函数计算

```bash
# 生成输入（自动添加真空层）
python3 setup_electronic.py --src CONTCAR --mode wf --vacuum 20 --ncore 8

# 运行 VASP
cd calc_electronic/wf_static
NP=16 EXE=vasp_std run_vasp.sh

# 后处理（计算功函数，绘制电势剖面）
python3 analyze_electronic.py --calcdir calc_electronic/wf_static --mode wf
```

**关键 INCAR 参数**:
| 参数 | 值 | 说明 |
|------|-----|------|
| LVHAR | .TRUE. | 输出 LOCPOT（静电势） |
| LDIPOL | .TRUE. | 启用偶极修正 |
| IDIPOL | 3 | z 方向偶极修正（slab） |
| ISMEAR | 0 | Gaussian 展宽 |

**输出**:
- `vacuum_potential_z.dat`: z 方向平面平均电势
- `wf_profile.png`: 电势剖面图
- 终端打印 E_F, V_vac, Φ

### DOS 计算（两步法）

```bash
# 生成输入（SCF + NSCF 两步）
python3 setup_electronic.py --src CONTCAR --mode dos --two_step

# 步骤 1: SCF 自洽
cd calc_electronic/dos_scf
NP=16 EXE=vasp_std run_vasp.sh

# 步骤 2: 拷贝 CHGCAR
cp CHGCAR ../dos_nscf/

# 步骤 3: NSCF 计算 DOS
cd ../dos_nscf
NP=16 EXE=vasp_std run_vasp.sh

# 后处理
python3 analyze_electronic.py --calcdir calc_electronic/dos_nscf --mode dos
```

**关键 INCAR 参数**:
| 参数 | 值 | 说明 |
|------|-----|------|
| ICHARG | 11 | 从 CHGCAR 读取电荷（NSCF） |
| LORBIT | 11 | 输出 PDOS（投影态密度） |
| NEDOS | 3000 | DOS 采样点数 |
| ISMEAR | -5 | 四面体法（半导体），金属用 0 |

**输出**:
- `dos_total.csv`: 能量-DOS 数据
- `dos.png`: DOS 图
- 推荐使用 `sumo` 或 `p4vasp` 进行更详细的 PDOS 分析

### setup_electronic.py 参数

```bash
python3 setup_electronic.py --help
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--src` | 必填 | 输入结构文件 |
| `--mode` | 必填 | wf 或 dos |
| `--outdir` | calc_electronic | 输出目录 |
| `--vacuum` | 20 | 真空层厚度 Å（wf 模式） |
| `--ncore` | - | NCORE 并行参数 |
| `--encut` | 500 | 截断能 eV |
| `--ediff` | 1e-6 | 收敛判据 |
| `--kpts_wf` | "8 8 1" | 功函数 K 点 |
| `--kpts_dos` | "12 12 1" | DOS K 点 |
| `--ismear_dos` | -5 | DOS 的 ISMEAR |
| `--two_step` | True | DOS 两步法 |

### 注意事项

- **功函数需要 slab 结构**：如果输入是 bulk，脚本会警告
- **VASP_PP_PATH**：需设置环境变量用于生成 POTCAR
  ```bash
  export VASP_PP_PATH=/path/to/potentials
  ```
- **DOS ISMEAR**: 半导体/绝缘体用 -5（四面体法），金属用 0

---

## 其他工具

### 检查计算状态

```bash
check_vasp.sh
```

### 清理文件

```bash
# 预览
DRYRUN=1 clean_vasp.sh

# 执行
clean_vasp.sh
```

### 配方换算

```bash
# 按目标原子数
python3 recipe_to_counts.py --target_atoms 200

# 按总质量
python3 recipe_to_counts.py --total_mass_g 1.0 --scale_to_atoms 5000
```

---

## 常见问题

**Q: oneAPI 加载失败？**
```bash
# 查看日志
cat /tmp/oneapi_setvars_$(id -u).log
```

**Q: 核数超限警告？**
```bash
# 自动下调（默认）
NP=32 run_vasp.sh

# 严格模式（报错退出）
STRICT_NP=1 NP=32 run_vasp.sh
```

**Q: 如何续算 AIMD？**
```bash
RESUME=1 NP=16 run_vasp.sh
```

**Q: 磁盘空间不足？**
```bash
# 强制继续
FORCE_DISK=1 run_vasp.sh

# 清理旧文件
clean_vasp.sh
```

**Q: ISYM 为什么强制为 0？**

AIMD 必须关闭对称性（ISYM=0），否则可能导致错误的轨迹或崩溃。脚本会强制覆盖 INCAR.base 中的 ISYM 设置。

**Q: MAXMIX 是什么？**

MAXMIX 控制电荷密度混合的历史长度，AIMD 中建议 40-80，有助于电子步收敛。

---

## 目录结构

```
~/vasp_scripts/
├── vasp_env.sh            # 环境配置（含自检）
├── run_vasp.sh            # 运行脚本（含续算/核数检查）
├── check_vasp.sh          # 状态检查
├── aimd_watch.sh          # AIMD 监控
├── aimd_msd.py            # MSD 分析
├── aimd_post.py           # 热力学后处理
├── clean_vasp.sh          # 文件清理
├── recipe.yaml            # 配方示例
├── recipe_validate.py     # 配方验证
├── recipe_to_counts.py    # 配方换算
├── make_incar_aimd.py     # AIMD INCAR 生成器
├── aimd_setup.sh          # AIMD 一键设置
├── setup_electronic.py    # 电子性质输入生成
├── analyze_electronic.py  # 电子性质后处理
├── setup_aimd_ase.py      # 大体系切割 AIMD cluster
└── README.md              # 本文档
```

## 版本信息

- VASP: 6.4.3
- Intel oneAPI: setvars.sh
- Python: 3.x (需要 numpy, pyyaml)

# VASP Scripts - 凝胶电解质计算工作流

针对 VASP (oneAPI + Intel MPI) 环境的计算辅助脚本集合，适用于凝胶电解质、AIMD 等计算任务。

---

## ⚠️ 使用边界（Scope）

> **在开始使用前，请理解 DFT-AIMD 的适用范围**

| 性质 | 推荐方法 | 说明 |
|------|----------|------|
| **输运性质（D, σ）** | 经典 MD (ns 级) | AIMD 时间尺度太短 (ps)，扩散系数仅供趋势参考 |
| **局域结构/配位** | DFT-AIMD | AIMD 优势：电子结构精确 |
| **反应趋势/过渡态** | DFT-AIMD / NEB | 反应能垒、过渡态搜索 |
| **电化学稳定窗口 (ESW)** | 片段氧化/还原能 | PBE 带隙不能直接当 ESW！ |
| **长程扩散 (ns 级)** | 经典 MD | AIMD 不适用 |

### 模型类型说明

- **bulk 模式**：周期性子胞，适合局域性质分析
- **cluster 模式**：有限真空簇，**不能用于 bulk 性质**，表面效应显著

---

## 🚨 常见陷阱（必读）

### 1. Cluster Trap（真空簇误用）
❌ **错误**：用 cluster 模式计算扩散系数并与 bulk 比较  
✅ **正确**：cluster 仅用于局域电子结构分析；扩散用 bulk 模式或经典 MD

### 2. False Diffusion Trap（伪扩散）
❌ **错误**：AIMD 几 ps 直接用 r(t)-r(0) 拟合 MSD 得到 D  
✅ **正确**：使用 `aimd_msd.py` v2.2（默认 MTO）：
  - Multiple Time Origins (MTO): 对每个 lag τ 平均所有起点
  - 检查 log-log 斜率 α ≈ 1（正常扩散）
  - 检查 D(t) 平台稳定
  - α < 0.8 表示亚扩散/caging，D 不可信

### 6. Bulk Density Trap（bulk 密度失控）
❌ **错误**：用 bounding box + buffer 定 bulk 子胞盒子 → 密度远低于凝胶  
✅ **正确**：`setup_aimd_ase.py` v2.2 按原体系密度反推 V = M_sub / ρ_orig

### 3. Langevin Gamma Trap（热浴干扰）
❌ **错误**：gamma=20 跑全程，扩散被抑制  
✅ **正确**：平衡段 gamma=10-20，生产段 gamma=1-5

### 4. Band Gap ≠ ESW Trap（带隙误用）
❌ **错误**：PBE 带隙 = 电化学稳定窗口  
✅ **正确**：ESW 需要片段氧化/还原能或反应自由能分析

### 5. Recipe Rounding Trap（凑整误差）
❌ **错误**：200 原子体系严格保持 wt%  
✅ **正确**：检查 counts_report.txt 误差；增大 target_atoms 减小误差

---

## 文件列表

| 文件 | 用途 |
|------|------|
| `vasp_env.sh` | VASP 运行环境配置（oneAPI、自检、线程控制） |
| `run_vasp.sh` | VASP 运行脚本（备份、日志、续算、核数检查） |
| `check_vasp.sh` | 检查计算状态（完成标志、能量、费米能级） |
| `aimd_watch.sh` | AIMD 实时监控（温度、离子步、能量） |
| `aimd_msd.py` | MSD 计算与扩散系数拟合 v2.2（MTO + α 判定 + 分段误差） |
| `aimd_post.py` | AIMD 热力学数据后处理（E0、T、F 导出 CSV） |
| `clean_vasp.sh` | 安全清理大文件（WAVECAR、CHGCAR 等） |
| `recipe.yaml` | 配方定义文件示例（8 类组分 + 模拟条件） |
| `recipe_validate.py` | 配方验证工具（校验 wt% 总和、温度、格式） |
| `recipe_to_counts.py` | 配方换算工具（wt% → 分子/原子数） |
| `make_incar_aimd.py` | AIMD INCAR 生成器（核心） |
| `aimd_setup.sh` | AIMD 一键设置脚本 |
| `setup_electronic.py` | 电子性质输入生成（功函数/DOS） |
| `analyze_electronic.py` | 电子性质后处理（功函数/DOS） |
| `setup_aimd_ase.py` | 从大体系切割 AIMD 子体系 v2.2（密度定盒 + 残基电荷） |
| `smoke_test.sh` | 功能验证脚本 |
| `examples/` | 示例文件目录 |

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

### 3. 功能验证

```bash
./smoke_test.sh
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

## 从大体系切割 AIMD 子体系（ASE）

### 概述

`setup_aimd_ase.py` 用于从大体系结构（GROMACS/Packmol 输出）中切割出一个可用于 VASP AIMD 的局部量子区域。

**v2.2 关键改进**：
- ✅ **按密度定盒子**：bulk 模式用 V = M_sub / ρ_orig（避免"低压气相"陷阱）
- ✅ **MIC 重成像**：选中原子按 MIC 位移重建坐标，保证空间连贯
- ✅ **残基电荷估计**：支持 TFSI（双三氟甲磺酰亚胺）等多原子离子
- ✅ **切断键检测**：警告可能的自由基/断链
- ✅ **键跳扩展**：`--bond_hops` 避免切断聚合物链
- ✅ **密度 meta**：输出 `density_original/target/achieved`

**物理原理**：
```
Bulk 盒子体积:
  V_target = M_sub / ρ_target
  
  其中 ρ_target 默认取原体系密度 ρ_orig
  这保证了子体系密度与原体系一致，物理合理

错误做法（旧版）:
  V = bbox + buffer  →  密度远低于实际，产生"低压气相"
```

**典型场景**：凝胶电解质体系有数千原子，AIMD 只能算几百原子，需要切割一个以目标离子为中心的小子体系。

### 基本用法

```bash
# bulk 模式（默认，周期性子胞，推荐用于凝胶电解质）
python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom Li --radius 8 --mode bulk --outdir aimd_bulk

# cluster 模式（真空簇，需显式指定）
python3 setup_aimd_ase.py --src equilibrated.pdb --center_atom Li --radius 8 --mode cluster --vacuum 20 --outdir aimd_cluster

# 保留完整分子 + 电荷中和
python3 setup_aimd_ase.py --src system.pdb --center_atom Li --radius 8 --selection molecule --neutralize nearest_counterions --outdir aimd_mol
```

> ⚠️ **重要**：默认为 bulk 模式。cluster 模式必须显式指定 `--mode cluster`。

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--src` | 必填 | 输入结构文件（.pdb/.gro/.xyz/.cif） |
| `--center_atom` | 必填 | 中心原子（索引或元素符号如 Li） |
| `--mode` | bulk | 模式: bulk（周期性）/ cluster（真空簇） |
| `--radius` | 8.0 | 切割半径 Å |
| `--selection` | sphere | 选择模式: sphere / molecule |
| `--bond_hops` | 0 | 键跳扩展步数（避免切断聚合物链） |
| `--density_g_cm3` | 原体系 | 目标密度 g/cm³（bulk 模式） |
| `--cell_shape` | scale_parent | 盒子形状: scale_parent / cubic |
| `--vacuum` | 20 | 真空层厚度 Å（仅 cluster 模式） |
| `--neutralize` | none | 电荷中和: none / nearest_counterions |
| `--charge_map_file` | - | 自定义残基电荷映射文件 |
| `--thermostat` | langevin | 恒温器: langevin / nose |
| `--gamma_1ps` | 10 | Langevin gamma (1/ps) |
| `--max_atoms` | 400 | 最大原子数限制 |
| `--temp` | 350 | AIMD 温度 K |
| `--steps` | 2000 | AIMD 步数 |
| `--potim` | 1.0 | 时间步长 fs（含 H 建议 0.5-1.0） |
| `--outdir` | aimd_sub | 输出目录 |

### 选择模式

- **sphere**：纯半径选择，最快
- **molecule**：半径内的原子所属分子整体保留，避免切断分子（需要结构文件包含分子信息，如 .pdb 的 residue）

> **建议**：凝胶电解质体系使用 `--selection molecule` 更物理合理

### 输出文件

```
aimd_Li8A/
├── POSCAR               # VASP 结构文件
├── INCAR                # AIMD 参数
├── KPOINTS              # K 点（Gamma-only）
├── POTCAR               # 赝势（如 VASP_PP_PATH 已设置）
├── cluster_visual.xyz   # 可视化文件（OVITO/VMD）
├── model_meta.json      # 元数据（含 density_*, 电荷, 警告）
├── selected_indices.txt # 选中的原子索引
└── cut_bonds_report.txt # 切断键报告（若有）
```

**model_meta.json 关键字段**:
```json
{
  "density": {
    "original_g_cm3": 1.18,
    "target_g_cm3": 1.18,
    "achieved_g_cm3": 1.17
  },
  "estimated_charge": 0,
  "cut_bonds": { "count": 0 },
  "warnings": []
}
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
python3 aimd_msd.py --specie Li --dt_fs 1.0 --t_skip_ps 2.0
```

### MSD 分析 (v2.2)

```bash
# 基本用法（默认 MTO）
python3 aimd_msd.py --specie Li --dt_fs 1.0

# 快速模式（stride=2 降低计算量）
python3 aimd_msd.py --specie Li --dt_fs 1.0 --stride 2

# 旧版兼容模式（单一时间原点，仅用于对比）
python3 aimd_msd.py --specie Li --dt_fs 1.0 --msd_method single_origin

# 指定拟合区间
python3 aimd_msd.py --specie Li --dt_fs 1.0 --t_skip_ps 2.0 --t_fit_start_ps 5.0

# 分段独立误差估计
python3 aimd_msd.py --specie Li --dt_fs 1.0 --n_blocks 4

# 关闭 unwrap 一致性检查（不推荐）
python3 aimd_msd.py --specie Li --dt_fs 1.0 --no_unwrap_check
```

**输出文件**：
- `msd_Li.dat`: MSD 数据（lag_ps, MSD_A2, n_samples）
- `D_running_Li.dat`: Running-D（D_ratio + D_deriv）
- `alpha_Li.dat`: log-log 斜率 α(t)
- `msd_report.txt`: 完整分析报告

**物理原理**：
```
MTO MSD 公式:
  MSD(τ) = ⟨|r(t₀+τ) - r(t₀)|²⟩_{t₀,ions}
  
  对每个 lag τ，平均所有时间起点 t₀ 和目标离子

α(t) 判定:
  α = d log(MSD) / d log(t)
  
  α ≈ 1: 正常扩散 (Fickian)
  α < 1: 亚扩散/受限 (caging, network constraint)
  α > 1: 超扩散/弹道 (early ballistic)
  
  只有 α ≈ 1 且 D(t) 平稳时，扩散系数才可信
```

**v2.2 关键改进**：

| 功能 | v2.1 | v2.2 |
|------|------|------|
| MTO | ✓ | ✓ + 智能 max_lag 默认 |
| 误差估计 | trajectory blocks | + 正确说明物理意义 |
| 计算效率 | - | **--stride 控制** |
| Unwrap | 跳变检测 | **+ |d|>0.5 一致性检查** |
| Vacuum 检测 | - | **自动识别 cluster** |

**新增/更新参数**：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--msd_method` | mto | mto / single_origin |
| `--stride` | 1 | MTO lag 步进（2/5 可降低计算量） |
| `--max_lag_ps` | min(T×0.5, 10) | 智能默认，避免噪声 |
| `--unwrap_check` | 启用 | 检测 |d|>0.5 分数坐标跳跃 |
| `--remove_com` | all | COM 漂移去除 |
| `--runningD` | both | ratio/derivative/both |
| `--alpha_window` | 21 | α(t) 滑窗大小 |
| `--block_mode` | trajectory_blocks | trajectory_blocks/bootstrap |
| `--seed` | - | Bootstrap 随机种子 |

**log-log 斜率 α(t) 解释**：

| α 值 | 状态 | 含义 |
|------|------|------|
| α ≈ 1 | diffusive | 正常扩散 |
| α < 0.8 | subdiffusive | 受限/亚扩散（caging） |
| α > 1.2 | superdiffusive | 弹道/超扩散（早期或漂移） |

> ⚠️ **重要**：只有 α ≈ 1 且 D(t) 稳定时，扩散系数才可信
> 
> ⚠️ 若检测到亚扩散（α < 0.8），说明 AIMD 时间内原子运动受限，扩散系数可能被高估

### make_incar_aimd.py 功能

- 从 recipe.yaml 读取模拟条件
- 温度自动转换：°C → K
- 强制设置 `ISYM = 0`（AIMD 必须）
- 设置 `MAXMIX`（默认 40）
- 添加 `LASPH = .TRUE.` 和 `ADDGRID = .TRUE.`
- 支持 INCAR.base 继承
- **v2.0**：自动检测含 H 体系并调整 POTIM
- **v2.0**：Langevin gamma 过大时警告
- **v2.0**：支持两段式 INCAR（平衡/生产）

### 两段式 AIMD（平衡/生产分离）

```bash
# 生成两段式 INCAR
python3 make_incar_aimd.py --two_stage

# 1. 平衡段（强控温）
cp INCAR.eq INCAR
NP=16 run_vasp.sh

# 2. 切换到生产段
cp CONTCAR POSCAR
cp INCAR.prod INCAR
NP=16 run_vasp.sh

# 3. 扩散分析只用生产段
python3 aimd_msd.py --specie Li --dt_fs 1.0
```

> ⚠️ **扩散系数只能用生产段数据**，平衡段 gamma 较大会抑制动力学

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

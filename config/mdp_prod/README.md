# 生产级 MDP 参数集 (Production MDP)

> 本目录包含用于**正式模拟/生产运行**的 GROMACS mdp 文件集合。

## 📋 设计理念

这套 mdp 文件专为**长时间生产模拟** (ns ~ 100 ns 级别) 设计，适用于：

- 聚合物电解质体系
- 凝胶电解质体系
- 锂电池/钠电池电解质模拟

与 `configs/mdp/` (smoke test) 的主要区别：

| 特性 | Smoke (`configs/mdp/`) | Production (`config/mdp_prod/`) |
|------|------------------------|----------------------------------|
| 目的 | 流程验证 | 正式统计采样 |
| 典型时长 | 100 ps ~ 1 ns | 10 ns ~ 100 ns |
| 收敛标准 | 宽松 | 严格 |
| 压力耦合 | Berendsen | Berendsen (预平衡) + P-R (生产) |
| 输出频率 | 可能过密 | 优化后的合理频率 |

## 📁 文件列表

| 文件 | 用途 | 默认时长 | 说明 |
|------|------|----------|------|
| `em.mdp` | 能量最小化 | ~50k 步 | 严格收敛 (emtol=500) |
| `nvt.mdp` | NVT 平衡 | 1 ns | 建立温度分布 |
| `npt_ber.mdp` | NPT 预平衡 | 2 ns | Berendsen 快速收敛密度 |
| `npt_pr.mdp` | NPT 严格平衡 | 5 ns | Parrinello-Rahman |
| `md_prod.mdp` | 生产 MD | 10 ns | 正式统计采样 |

## 🔄 推荐流程

```
em.mdp → nvt.mdp → npt_ber.mdp → npt_pr.mdp → md_prod.mdp
  │         │           │            │            │
能量最小化  温度平衡    密度快速收敛   严格NPT平衡   正式采样
```

### 流程说明

1. **em.mdp**: 消除初始结构的高能构象和原子重叠
2. **nvt.mdp**: 在恒定体积下建立合理的温度分布
3. **npt_ber.mdp**: 使用 Berendsen 压控快速收敛密度 (⚠️ 不用于统计)
4. **npt_pr.mdp**: 使用 Parrinello-Rahman 进行严格 NPT 平衡
5. **md_prod.mdp**: 生产段，用于计算统计性质

### ⚠️ 重要提示

- **Berendsen 压控不产生正确的 NPT 系综**，仅用于预平衡
- 正式统计必须使用 **Parrinello-Rahman** (npt_pr 或 md_prod)
- 直接从 em → md_prod 跳过平衡步骤会导致体系不稳定

## ⏱️ 修改模拟时长

所有 mdp 文件使用 `dt = 0.002 ps` (2 fs) 时间步。修改时长只需改 `nsteps`:

```
时长 = nsteps × dt
```

### 换算表

| 目标时长 | nsteps |
|----------|--------|
| 1 ns | 500,000 |
| 2 ns | 1,000,000 |
| 5 ns | 2,500,000 |
| 10 ns | 5,000,000 |
| 20 ns | 10,000,000 |
| 50 ns | 25,000,000 |
| 100 ns | 50,000,000 |

### 示例：将 md_prod 从 10 ns 改为 50 ns

```bash
# 方法 1: 直接编辑 md_prod.mdp
sed -i 's/nsteps.*=.*5000000/nsteps = 25000000/' config/mdp_prod/md_prod.mdp

# 方法 2: 使用 grompp 的 -maxh 参数控制运行时间
gmx grompp -f md_prod.mdp -c npt_pr.gro -t npt_pr.cpt -p topol.top -o md_prod.tpr
gmx mdrun -deffnm md_prod -maxh 72  # 最多运行 72 小时
```

## 🌡️ 修改温度/压力

温度和压力参数在每个 mdp 文件中都有明确标注：

```mdp
; 温度
ref_t = 298           ; 修改目标温度 (K)
gen_temp = 298        ; nvt.mdp 中的初始速度温度

; 压力
ref_p = 1.0           ; 修改目标压力 (bar)
compressibility = 4.5e-5  ; 等温压缩系数
```

### 不同体系的压缩系数建议

| 体系类型 | compressibility (bar⁻¹) |
|----------|-------------------------|
| 水溶液 | 4.5e-5 |
| 有机溶剂 | 4.5e-5 ~ 6e-5 |
| 聚合物溶液 | 3e-5 ~ 4.5e-5 |
| 聚合物固体 | 1e-5 ~ 2e-5 |

## 📊 输出文件说明

### 默认输出频率

| 参数 | 值 | 实际频率 | 用途 |
|------|-----|---------|------|
| `nstxout-compressed` | 5000 | 10 ps | 轨迹 (xtc) |
| `nstenergy` | 500 | 1 ps | 能量 (edr) |
| `nstlog` | 500 | 1 ps | 日志 (log) |

### 文件大小估算 (10,000 原子体系)

| 文件类型 | 每 10 ns 大小 |
|----------|---------------|
| .xtc | ~100 MB |
| .edr | ~50 MB |
| .log | ~10 MB |

### 需要更精细轨迹？

如需更精细的 MSD 分析，可修改 `nstxout-compressed`:

```mdp
nstxout-compressed = 2500    ; 5 ps/帧 (文件大小翻倍)
nstxout-compressed = 1000    ; 2 ps/帧 (文件大小 5x)
```

## 🔧 参数一致性

所有 mdp 文件保持以下参数一致，确保平滑过渡：

```mdp
; 非键相互作用
cutoff-scheme = Verlet
rcoulomb = 1.2
rvdw = 1.2
coulombtype = PME
vdwtype = Cut-off
vdw-modifier = Potential-shift
DispCorr = EnerPres

; 约束
constraints = h-bonds
constraint-algorithm = LINCS
lincs-order = 4

; 周期性
pbc = xyz
```

## 🚀 快速开始

### 使用 run_workflow.sh (推荐)

```bash
# 设置使用生产 mdp
export MDP_DIR=config/mdp_prod

# 运行完整流程
./scripts/run_workflow.sh
```

### 手动运行

```bash
MDP=config/mdp_prod

# 1. 能量最小化
gmx grompp -f $MDP/em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

# 2. NVT 平衡
gmx grompp -f $MDP/nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# 3. NPT 预平衡 (Berendsen)
gmx grompp -f $MDP/npt_ber.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt_ber.tpr
gmx mdrun -deffnm npt_ber

# 4. NPT 严格平衡 (P-R)
gmx grompp -f $MDP/npt_pr.mdp -c npt_ber.gro -t npt_ber.cpt -p topol.top -o npt_pr.tpr
gmx mdrun -deffnm npt_pr

# 5. 生产 MD
gmx grompp -f $MDP/md_prod.mdp -c npt_pr.gro -t npt_pr.cpt -p topol.top -o md_prod.tpr
gmx mdrun -deffnm md_prod
```

## 📚 参考

- GROMACS 2023 Manual: https://manual.gromacs.org/2023/
- MDP options: https://manual.gromacs.org/current/user-guide/mdp-options.html

---

*生成于 gel_packmol 项目*


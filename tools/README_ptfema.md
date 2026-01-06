# PTFEMA 链构建工具

构建 PTFEMA（聚三氟乙基甲基丙烯酸酯 / Poly(2,2,2-trifluoroethyl methacrylate)）线性链的工具。

## 化学说明

- **单体**: TFEMA (2,2,2-三氟乙基甲基丙烯酸酯)
- **聚合位点**: C=C 双键（甲基丙烯酸酯类聚合）
- **重复单元**: –CH₂–C(CH₃)(COOCH₂CF₃)–
- **主链结构**: 碳-碳单键主链

## 安装依赖

```bash
conda install -c conda-forge rdkit
pip install python-docx  # 可选，用于解析 docx 坐标
```

## 使用方法

### 基本命令

```bash
cd ~/gel_packmol
python tools/build_ptfema.py --monomer mol2/TFEMA.mol2 --n 10 --cap H
```

### 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--monomer` | 单体 MOL2 文件路径 | (必需) |
| `--n` | 聚合度 | 10 |
| `--cap` | 端基类型: `H` 或 `Me` | H |
| `--outprefix` | 输出文件前缀 | PTFEMA |
| `--outdir` | 输出目录 | 与单体同目录 |
| `--docx` | 补充材料 docx 文件路径 | (可选) |

## 常用命令示例

### 1. H 端基封端（最常用）

```bash
python tools/build_ptfema.py \
    --monomer mol2/TFEMA.mol2 \
    --n 10 \
    --cap H \
    --outprefix PTFEMA
```

输出:
- `mol2/PTFEMA_repeat_capped.mol2` - 重复单元 + 端基封端
- `mol2/PTFEMA_10.mol2` - 10 聚体 MOL2
- `mol2/PTFEMA_10.pdb` - 10 聚体 PDB

### 2. 甲基端基封端

```bash
python tools/build_ptfema.py \
    --monomer mol2/TFEMA.mol2 \
    --n 10 \
    --cap Me \
    --outprefix PTFEMA_Me
```

输出:
- `mol2/PTFEMA_Me_repeat_capped.mol2`
- `mol2/PTFEMA_Me_10.mol2`
- `mol2/PTFEMA_Me_10.pdb`

### 3. 使用 DOCX 坐标模板

```bash
python tools/build_ptfema.py \
    --monomer mol2/TFEMA.mol2 \
    --n 10 \
    --cap H \
    --docx "supplementary data 1.docx" \
    --outprefix PTFEMA
```

输出（额外）:
- `mol2/PTFEMA_from_docx.xyz` - 从 docx 提取的坐标
- `mol2/PTFEMA_from_docx.mol2` - 推断键后的 MOL2（如果成功）
- `mol2/PTFEMA_from_docx.pdb` - 推断键后的 PDB（如果成功）

## 输出文件说明

### PTFEMA_repeat_capped.mol2

重复单元 + 端基封端的小模型，用于：
- 力场参数化
- 电荷拟合
- 基准计算

### PTFEMA_N.mol2 / PTFEMA_N.pdb

N 聚合度的完整链，用于：
- MD 模拟
- 构象分析
- 可视化

## 验证输出

脚本运行时会打印：

1. **聚合双键识别**: CH2 端碳和取代端碳的原子索引
2. **片段数**: 必须 = 1（单一连通分子）
3. **端基类型**: H 或 Me
4. **主链键长**: min/max/avg（应在 ~1.50-1.60 Å）
5. **几何检查**: 最小非键距离（应 > 0.8 Å）

## 故障排除

### "无法识别聚合双键"

确保输入的 MOL2 文件：
- 包含 `@<TRIPOS>BOND` 部分
- C=C 双键的键阶标记为 `2`
- 是甲基丙烯酸酯类结构

### "ETKDG 失败"

分子嵌入困难时会自动使用随机坐标。如果最终结构有问题，可以：
1. 手动用可视化软件调整
2. 用外部工具（如 Avogadro）优化

### docx 解析失败

需要安装 `python-docx`:

```bash
pip install python-docx
```

## 技术细节

### 聚合双键识别算法

1. 遍历所有键阶=2 的键
2. 筛选 C-C 双键（排除 C=O 羰基）
3. 判断哪个碳连着羰基碳（取代端）
4. 另一个碳为 CH2 端

### 链构建策略

1. 构建 SMILES 表示的链
2. RDKit 解析 SMILES
3. ETKDG 嵌入 3D 坐标
4. MMFF/UFF 力场优化

### 文件格式

- MOL2: 自定义读写器，不依赖 OpenBabel
- PDB: 使用 RDKit 的 `MolToPDBFile`
- XYZ: 简单文本格式


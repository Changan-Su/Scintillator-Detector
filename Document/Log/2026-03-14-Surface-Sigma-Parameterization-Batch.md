# Log: Surface Sigma 参数化与 run_batch 集成

**日期**: 2026-03-14  
**项目**: Scintillator-Detector-Continuous  

---

## 开发内容概览

1. **Surface Sigma 宏可调**：将晶体光学表面粗糙度 `Surface_Sigma` 从 `DetectorConstruction.cc` 硬编码改为通过 `/detector/surfaceSigma` 宏命令在运行时设定。
2. **run_batch Sigma 循环**：在 `run_batch.bat` 与 `run_batch.sh` 中增加 Sigma 参数扫描支持，可全自动按设定数值循环跑结果。
3. **文档更新**：`Document/Notes/Note-2026-03-14-Surface-Sigma-Parameterization.md` 记录实现方案与知识点；`BATCH_USAGE.md` 增加 Sigma 扫描示例。

---

## 实现方式

### 1) DetectorConstruction / DetectorMessenger

- **`include/DetectorConstruction.hh`**：新增 `fSurfaceSigma`（默认 0.5）、`SetSurfaceSigma(G4double)`、`GetSurfaceSigma()`。注意使用 `G4double`（小写 d），`G4Double` 未定义会导致编译失败。
- **`include/DetectorMessenger.hh`**：新增 `fSurfaceSigmaCmd` 指针声明。
- **`src/DetectorMessenger.cc`**：创建 `/detector/surfaceSigma` 命令（`G4UIcmdWithADouble`，范围 0～1），析构中 `delete`，`SetNewValue` 中调用 `SetSurfaceSigma`。
- **`src/DetectorConstruction.cc`**：`Construct()` 中 `Surface_Sigma = fSurfaceSigma`，用于 `crystalsurface->SetSigmaAlpha(Surface_Sigma)`。

### 2) run_batch 集成

- **配置项**：`LOOP_SURFACE_SIGMA`、`SIGMA_START`/`SIGMA_END`/`SIGMA_STEP`、`DEFAULT_SURFACE_SIGMA`、`NAME_INCLUDE_SURFACE_SIGMA`。
- **generate_geometry_mac**：增加第 11 个参数，输出 `/detector/surfaceSigma <值>`。
- **嵌套循环**：Sigma 作为最外层循环，与 Nx/Ny/Nz/Gap/Size 等组合。
- **文件夹命名**：`_Sigma0p7` 表示 sigma=0.7（小数点用 `p` 替代）。

---

## 涉及文件

| 文件 | 修改内容 |
|------|----------|
| `include/DetectorConstruction.hh` | `fSurfaceSigma`、setter/getter |
| `include/DetectorMessenger.hh` | `fSurfaceSigmaCmd` 声明 |
| `src/DetectorMessenger.cc` | 命令创建、析构、SetNewValue |
| `src/DetectorConstruction.cc` | 使用 `fSurfaceSigma` |
| `run_batch.bat` | Sigma 配置、循环、generate_geometry_mac |
| `run_batch.sh` | 同上 |
| `Document/Notes/Note-2026-03-14-Surface-Sigma-Parameterization.md` | 实现与知识点笔记 |
| `BATCH_USAGE.md` | 示例 4：Surface Sigma 扫描 |

---

## 使用方式

**单次运行**：在 `geometry.mac` 中增加 `/detector/surfaceSigma 0.7`。

**批量扫描**：设置 `LOOP_SURFACE_SIGMA=true`，配置 `SIGMA_START`/`SIGMA_END`/`SIGMA_STEP`，运行 `run_batch.bat` 或 `run_batch.sh`。

---

## 增量更新：Histo10 结果 CSV 落盘（重建坐标）

### 本次开发内容

- 在 `Histo10_Cubic.py` 中，除了原有 `merged_event.csv`、`merged_face_jk.csv` 和热力图外，新增输出 `reconstructed_position.csv`。
- 新增 CSV 用于保存当前配置目录的重建结果，便于后续批量对比与二次分析，不需要再从图片标题/图下注释手工抄录。

### 实现方式

- 新增函数 `compute_face_totals(agg)`：从 `Face,j,k,Count` 聚合表中计算六个面的总光子数（Face 0..5）。
- 在 `process_one_config(...)` 中：
  - 先基于 `agg` 计算六面总计数；
  - 调用 `compute_reconstructed_position(...)` 计算归一化重建坐标；
  - 调用 `compute_reconstructed_position_2(...)` 计算 sqrt 版本重建坐标；
  - 将以上结果写入 `reconstructed_position.csv`。
- 输出 CSV 字段包含：
  - 配置名：`ConfigName`
  - 六面总计数：`N_plusX/N_minusX/N_plusY/N_minusY/N_plusZ/N_minusZ`
  - 重建坐标：`x_rec/y_rec/z_rec`
  - sqrt 公式重建坐标：`x_rec_sqrt/y_rec_sqrt/z_rec_sqrt`
  - 总光子数：`PhotonTotal`

### 说明

- 这次改动是新增结果文件，不会改变原始合并数据和热力图生成逻辑。
- 先前热力图底部仅显示重建坐标文字，本次将同类结果结构化保存为 CSV，可直接用于后续统计脚本读取。

---

## 增量更新：Results 目录命名简化 + metadata.csv 宏参数落盘

### 本次开发内容

- **Results 子目录命名简化**：从原来超长格式（`_Nx_Ny_Nz_Gap_Size_...`）改为纯时间戳（`20260314_170042_581`）。
- **`/results/prefix` 宏命令**：新增宏命令，可在宏文件中设置输出目录前缀。前缀非空时格式为 `<prefix>_<timestamp>`；前缀为空时退化为纯时间戳。
- **`metadata.csv`**：每个 Results 子目录新增 `metadata.csv`，以 `Key,Value` 两列格式记录本次仿真的宏配置参数，完全替代目录名中丢失的几何信息。

### 实现方式

- **`include/DetectorConstruction.hh`**：新增 `fResultsPrefix`（`G4String`，默认空）、`GetResultsPrefix()`、`SetResultsPrefix()`。
- **`include/DetectorMessenger.hh` / `src/DetectorMessenger.cc`**：注册 `/results/prefix` 命令（`G4UIcmdWithAString`），`SetNewValue` 中调用 `SetResultsPrefix`。
- **`include/PrimaryGeneratorAction.hh`**：增加 `GetSourceMode()`、`GetSourceDistribution()`、`GetFpSource()` 三个只读 getter，供 metadata 写出使用。
- **`src/HistoManager.cc`**：
  - 目录名构建改为 `prefix + "_" + timestamp`（prefix 为空则只 timestamp）。
  - 首次创建目录时（master 线程）写出 detector 参数行。
  - 第一个 worker 线程拿到 `PrimaryGeneratorAction` 实例后追加写入 source 参数行（解决 MT 模式下 master 线程无法访问 `PrimaryGeneratorAction` 的问题）。
  - 全局 flag `gMetadataSourceWritten` 保证 source 参数只写一次。

### metadata.csv 字段清单

| Key | 来源 |
|-----|------|
| `detector.surfaceSigma` | `/detector/surfaceSigma` |
| `detector.arrayNx/Ny/Nz` | `/detector/arrayNx/Ny/Nz` |
| `detector.crystalGap_mm` | `/detector/crystalGap` |
| `detector.crystalSize_mm` | `/detector/crystalSize` |
| `detector.crystalSizeY_mm` | `/detector/crystalSizeY` |
| `detector.fillterRatioY/Z` | `/detector/fillterRatioY/Z` |
| `detector.fillterPosRatioY/Z` | `/detector/fillterPosRatioY/Z` |
| `results.prefix` | `/results/prefix` |
| `source.mode` | `/source/mode` |
| `source.distribution` | `/source/distribution` |
| `source.fp_source_x/y/z_mm` | `/source/fp_source` |

### 使用方式

在宏文件中加入（`/run/initialize` 之前）：
```
/results/prefix MyTag
```
生成目录：`Results/MyTag_20260314_170057_934/`；不写则：`Results/20260314_170042_581/`。

---

## 增量更新：Histo10 输出目录同步复制 metadata.csv

### 本次开发内容

- 在 `Histo10_Cubic.py` 的单配置处理流程中，新增从 `Results/<配置目录>/metadata.csv` 到 `Output/<批次>/<配置目录>/metadata.csv` 的自动复制。
- 目的：让后处理输出目录中的图表与汇总 CSV 直接携带同目录参数元数据，避免分析时在 `Results` 与 `Output` 间来回查找。

### 实现方式

- 在 `process_one_config(...)` 开头增加：
  - `metadata_src = config_dir / "metadata.csv"`
  - 若源文件存在，则 `shutil.copy2(metadata_src, metadata_dst)` 复制到输出目录。
- 在处理完成日志中新增 `metadata.csv` 路径打印（仅复制成功时显示）。

### 影响说明

- 不影响既有 `merged_event.csv`、`merged_face_jk.csv`、`reconstructed_position.csv`、`sipm_6faces_heatmap.png` 的生成逻辑。
- 当某个配置目录没有 `metadata.csv` 时，脚本保持兼容，继续正常处理该目录。

---

## 增量更新：新增仅扫 source position + surfaceSigma 的批处理脚本

### 本次开发内容

- 新增脚本 `run_batch_sigma_position.bat`，用于只循环两类参数：
  - `/detector/surfaceSigma`
  - `/source/fp_source x y z`
- 其它参数保持固定（通过 `geometry.mac` 与脚本固定配置注入），避免沿用通用大脚本时需要关闭大量开关。

### 实现方式

- 在每组参数组合下动态生成临时宏 `run_tmp_sigma_position.mac`，内容包括：
  - `/control/execute geometry.mac`
  - `/detector/surfaceSigma <sigma>`
  - `/results/prefix <自动前缀>`
  - `/source/mode`、`/source/distribution`、`/source/fp_source x y z`
  - `/run/beamOn <N>`
- 前缀自动编码规则：`S<sigma>_X<x>_Y<y>_Z<z>`（小数点转 `p`、负号转 `m`），便于快速识别结果目录参数来源。
- 支持 `DRY_RUN=true` 干跑模式，仅打印组合与前缀，不执行 Geant4。

### 说明

- 脚本已完成干跑验证，确认负数坐标（如 `y=-7.5`）和小数步进可正确生成参数列表。
- 该脚本与现有 `run_batch.bat` 并存：前者用于“只扫位置+sigma”的高频实验，后者用于全参数组合扫描。

---

## 增量更新：新增 Python 版 position+sigma 扫描脚本

### 本次开发内容

- 新增 `run_batch_sigma_position.py`，功能与 `run_batch_sigma_position.bat` 对齐，但可读性更高、参数传递更直观。
- 支持仅循环：
  - `--sigma-start/--sigma-end/--sigma-step`
  - `--x-start/--x-end/--x-step`
  - `--y-start/--y-end/--y-step`
  - `--z-start/--z-end/--z-step`
- 每组参数自动生成临时宏并设置：
  - `/detector/surfaceSigma`
  - `/results/prefix`
  - `/source/fp_source x y z`
  - `/run/beamOn`

### 使用方式（示例）

```bash
python run_batch_sigma_position.py --dry-run
python run_batch_sigma_position.py --sigma-start 0.3 --sigma-end 0.7 --sigma-step 0.1 --x-start 8 --x-end 12 --x-step 2
```

### 说明

- `--dry-run` 仅打印参数组合和前缀，不执行 Geant4。
- 前缀自动编码格式：`S..._X..._Y..._Z...`（小数点转 `p`，负号转 `m`）。

---

## 增量更新：Histo10 同时输出多种可行重建算法

### 本次开发内容

- 将 `Histo10_Cubic.py` 中原先的两种重建结果扩展为一组统一输出的多算法框架。
- `reconstructed_position.csv` 现在会对每个配置目录输出多行，每行对应一种算法：
  - `linear_scaled`
  - `sqrt_scaled`
  - `log_scaled`
  - `centroid_mean`
  - `centroid_weighted`
  - `hybrid_sqrt_centroid`

### 实现方式

- 新增 `read_half_lengths_cm(...)`：优先从 `metadata.csv` 读取晶体尺寸，统一换算为 cm 半长。
- 新增三类“只用六面总光子数”的算法：
  - `linear_scaled`
  - `sqrt_scaled`
  - `log_scaled`
- 新增两类“利用每个面 4×4 分布质心”的算法：
  - `centroid_mean`
  - `centroid_weighted`
- 新增混合算法：
  - `hybrid_sqrt_centroid`（`sqrt_scaled` 与 `centroid_weighted` 各占 50%）
- 通过 `compute_reconstruction_rows(...)` 统一收集全部算法结果，后续若继续加算法，只需补函数并加入此列表。

### 热力图显示

- 热力图底部不再固定显示旧的单一公式结果，改为显示一个参考算法：
  - `hybrid_sqrt_centroid`
- 其余完整算法结果保存在 `reconstructed_position.csv` 中，供 `analyze_position_accuracy.py` 自动读取和评估。

### 增量更新（可读性）：差分三类算法改为“整函数展开”

- **原因**：原先 `linear_scaled` / `sqrt_scaled` / `log_scaled` 共用一个 `axis_component_from_pair(..., transform)`，公式上正确但阅读时需跳转到通用函数才能看清每轴在算什么。
- **改动**：删除该通用函数；在 `reconstruct_linear_scaled`、`reconstruct_sqrt_scaled`、`reconstruct_log_scaled` 三个函数内分别完整写出 +X/-X、+Y/-Y、+Z/-Z 的变量名与分母判零逻辑，并在 docstring 中写出对应数学式。**数值行为与改前一致**（仍为 `half * (F(p)-F(n))/(F(p)+F(n))`，无信号轴为 0）。

### 增量更新：half_side_ratio（按半空间选取 SiPM 格子求和）

- **内容**：新增算法 `half_side_ratio`：对每一轴用「+ 侧 / − 侧」两组 **(Face, j, k)** 列表（`HALF_SIDE_X_PLUS_CELLS` 等）对 `sum_sipm_cells` 求和，再按 `half * (N+−N−)/(N++N−)` 映射到 cm（与整面对比的 `linear_scaled` 同形，但分子分母来自选定格子而非整面总和）。
- **默认 X 侧**：+X 面上一行 4 格 + 四个侧面靠 ±X 各 2 格（共 8 格）；Y/Z 为对称默认，可只改列表不改比值函数。
- **辅助**：`sum_sipm_cells(matrices, cells)` 供后续自定义半空间定义时复用。

### 增量更新：half_side_ratio 半侧定义改为“主面整面 + 侧面半面”

- 根据后续确认，`half_side_ratio` 的默认半侧定义已调整为：某轴正侧 = **正向主面整面 4x4** + 四个相关侧面中**靠该正侧的 2x4 半面**；负侧同理。
- 例如 X 轴：`+X` 侧现在是 `Face 0` 的全部 16 格，再加 `Face 2/3/4/5` 上靠 `+X` 的 `j=2,3` 两列（每面 8 格）；`-X` 侧对应 `Face 1` 全部 16 格 + 四个侧面的 `j=0,1` 两列。
- Y/Z 轴按各自局部坐标的同样规则对称推广；比值公式本身不变，仍为 `half * (N+−N−)/(N++N−)`。

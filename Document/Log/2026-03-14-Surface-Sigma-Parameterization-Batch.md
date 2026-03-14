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

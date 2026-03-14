# Log: GAGG 唯一闪烁体、SiPM 六面 4×4、几何参数开关与可视化修复

**日期**: 2026-02-21  
**项目**: Scintillator-Detector-Continuous  

---

## 开发内容概览

**几何与可视化（2026-02-21）**  
1. **闪烁体仅用 GAGG**：移除 YSO，晶体与 Fillter 均使用 GAGG 材料及 GAGG-HL 光学/闪烁属性。  
2. **晶体阵列参数来源开关**：新增 bool 开关，可选“geometry.mac/UI 参数”或“手动参数”，并配套 `/detector/manual*` 等命令。  
3. **SiPM 六面 4×4 排布**：每个晶体 6 个面各贴 4×4 SiPM（像素 6 mm，间隙 0.2 mm），共 6×16=96 个 SiPM/晶体。  
4. **Qt 可视化六面显示修复**：通过唯一物理体积名与旋转矩阵常驻，使 6 个面的 SiPM 在 Qt/OpenGL 中均正确显示。  
5. **编译与编码**：关键注释改为英文、显式类型与括号匹配修正，消除 MSVC 代码页 936 下的解析错误。

**6 面 SiPM 统计与输出（2026-02-21～22）**  
6. **6 面 SiPM 统计**：按 SiPMBlockID（crystalId×6×16 + face×16 + j×4 + k）统计光子，输出事件明细 CSV（PhotonFaceBlockEvent）；SteppingAction/EventAction/HistoManager 联动。  
7. **旧统计 CSV 停用**：PhotonLRPerRod、PhotonLeft/Right、PhotonDepthNtuple 等创建与填充已注释化，仅保留 PhotonFaceBlockEvent。  
8. **BlockTotal 移除**：运行总计改为由 Python 脚本合并线程 CSV 后计算；RunAction::EndOfRunAction 中补充 Save() 保证 ntuple 落盘。  
9. **热力图脚本**：根目录 `plot_sipm_6faces_heatmap.py` 支持单 CSV 或 `--dir` 合并线程；`Histo10_Cubic.py` 实现 Histo9 风格批处理（遍历 Results 下配置、输出合并 CSV + 聚合 CSV + 6 面热力图到 `Output/SiPM6_Output_<时间戳>/<配置名>/`）。  
10. **输出路径说明**：Geant4 CSV 写于 **CWD + Results/<时间戳>_.../**（如从 build 启动则位于 build/Results）；Python 分析输出在 **CWD 下 Output/SiPM6_Output_<时间戳>/**，可用 `--output` 指定。

---

## 实现方式

### 1) GAGG 唯一闪烁体

- **文件**: `src/DetectorConstruction.cc`
- 保留 GAGG 材料定义（Gd3Al2Ga3O12，密度 6.6 g/cm³）及 GAGG-HL 材料属性表（折射率 1.91，光产额 54 000 ph/MeV，衰减 150 ns，发射峰约 530 nm）。
- 删除原 YSO 用 MPT 块及 `YSO->SetMaterialPropertiesTable(...)`。
- `logicCrystal`、`logicFillter` 的 `G4LogicalVolume` 构造均改为使用 `GAGG` 材料。
- 晶体与 Fillter 的闪烁与光学行为统一由 GAGG 的 `gagg_mt` 决定。

### 2) 晶体阵列参数来源开关（geometry.mac vs 手动）

- **文件**: `include/DetectorConstruction.hh`
  - 新增 `fUseGeometryMac`（默认 true）及一套手动参数：`fManual_nx/ny/nz`、`fManual_crystal_gap`、`fManual_crystal_l/ly`、`fManual_fillter_ratio_*`、`fManual_fillter_pos_*`。
  - 新增对应 setter：`SetUseGeometryMac`、`SetManualArrayNx/Ny/Nz`、`SetManualCrystalGap/Size/SizeY`、`SetManualFillterRatio*`、`SetManualFillterPosRatio*`。
- **文件**: `include/DetectorMessenger.hh`、`src/DetectorMessenger.cc`
  - 新增命令：`/detector/useGeometryMac`（bool）、`/detector/manualArrayNx/Ny/Nz`、`/detector/manualCrystalGap`、`/detector/manualCrystalSize/SizeY`、`/detector/manualFillterRatioY/Z`、`/detector/manualFillterPosRatioY/Z`。
- **文件**: `src/DetectorConstruction.cc`
  - 在晶体阵列构建前，根据 `fUseGeometryMac` 选择使用“geometry.mac 参数”（`fPar_*`、`fCrystal_gap`、`fcrystal_l/ly`、Fillter 比例）或“手动参数”（`fManual_*`），再计算 `Crystal_nx/ny/nz`、`Crystal_gap`、`crystal_l/ly` 及 Fillter 比例，后续几何与 SiPM 放置均基于此组变量。

### 3) SiPM 六面 4×4 排布

- **文件**: `src/DetectorConstruction.cc`
- **单像素与间隙**：`sipm_l = 6 mm`，`SiPm_gap = 0.2 mm`，每面 4×4 像素，中心距 `pitch = 6.2 mm`（对应 25 mm 晶体单边 4 像素 + 5 间隙）。
- **几何**：SiPM 为薄盒，厚度 `sipm_t = 0.6 mm` 沿法向贴晶体表面，6×6 mm 在面内；`G4Box` 半边长顺序为 (厚度/2, 边长/2, 边长/2)，放置时通过旋转使厚度方向对准各面法向。
- **6 个面**：face 0~5 对应 +X、-X、+Y、-Y、+Z、-Z；每个面中心在晶体外表面外偏 `sipm_t/2`，面内用 `(j-1.5)*pitch`、`(k-1.5)*pitch` 做 4×4 网格偏移。
- **旋转**：+X 不旋转（nullptr）；-X 绕 Y 转 180°；±Y 绕 Z 转 ±90°；±Z 绕 Y 转 ±90°。旋转矩阵在循环外创建 5 个（`rotSiPM_mX/pY/mY/pZ/mZ`）并复用，**不在循环内 delete**，避免几何树中旋转失效导致仅一面显示。
- **Copy 号与命名**：`copyNoSiPM = crystalId*6*16 + face*16 + j*4 + k`；物理体积名 `SiPM_%d`（copyNoSiPM）保证每个 SiPM 唯一名，便于 Qt 场景树与可视化正确显示 6 面。

### 4) 编译与编码修复

- **文件**: `src/DetectorConstruction.cc`
- 将 GAGG、晶体阵列、Fillter、SiPM 及三重循环与 6 面 SiPM 段落中的**中文注释改为英文**，避免在 MSVC 代码页 936 下多字节字符被误解析为 `}` 等，导致“未声明的标识符”“非法 else”“缺少类型说明符”等连锁错误。
- `logicCrystal`、`solidCrystal` 改为显式类型 `G4LogicalVolume*`、`G4Box*`，避免在 YSO 移除后类型推断异常。
- 保持 `Construct()` 内大括号匹配正确，`fScoringVolume`、`fCrystal_*`、`flogicSiPM` 赋值与 `return physWorld` 均在函数内。

---

## 验证与注意事项

- **构建**：`cmake --build build --config Release` 通过，生成 `exampleB1.exe`；若出现 C4819，可将源文件存为“带 BOM 的 UTF-8”或保持当前英文注释。
- **运行**：geometry.mac 或 UI 设置晶体参数后，`/detector/update` 重建几何；手动模式时先 `/detector/useGeometryMac false`，再设 `/detector/manual*` 后 update。
- **可视化**：Qt 中 6 个面的 SiPM 应均可见；若曾出现“仅一面 16 个”，系旋转矩阵在循环内被 delete 导致，已通过“循环外创建、复用、不 delete”修复。

---

## 对既有开发内容的修正说明

- 此前 SiPM 仅在 ±X 两列放置且依赖已注释的 `sipm_l_ratio`，现改为每晶体 6 面 × 4×4，并统一 6 mm / 0.2 mm 规格。
- 晶体阵列参数曾仅来自 geometry.mac；现可通过 `fUseGeometryMac` 与 `/detector/manual*` 使用一套独立手动参数，便于脚本或固定配置测试。

---

## 涉及文件一览

| 文件 | 修改要点 |
|------|----------|
| `src/DetectorConstruction.cc` | GAGG 唯一、几何分支、SiPM 6 面 4×4、旋转常驻、注释/类型/括号 |
| `include/DetectorConstruction.hh` | `fUseGeometryMac`、`fManual_*` 及 setter |
| `src/DetectorMessenger.cc` | `/detector/useGeometryMac`、`/detector/manual*` 命令 |
| `include/DetectorMessenger.hh` | 上述命令成员指针 |

---

## 增量更新：6面 SiPM 统计与热力图（2026-02-21）

### 本阶段开发内容

1. **统计口径从 Left/Right 扩展到 6 面 SiPM 块**
   - 以 `copyNoSiPM = crystalId*6*16 + face*16 + j*4 + k` 作为唯一 `SiPMBlockID`。
   - 命中 SiPM 时按 `SiPMBlockID` 计数，不再使用旧版 `endTag(Left/Right)` 逻辑。

2. **新增两类 CSV 输出**
   - **事件明细**：`PhotonFaceBlockEvent`  
     列：`EventID, CrystalID, iy, iz, Face, j, k, SiPMBlockID, PhotonCount`
   - **运行总计**：`PhotonFaceBlockTotal`  
     列：`CrystalID, iy, iz, Face, j, k, SiPMBlockID, PhotonCountTotal`

3. **新增根目录热力图脚本**
   - 文件：`plot_sipm_6faces_heatmap.py`
   - 功能：输入事件明细或运行总计 CSV，聚合为 6 个面的 `4x4` 计数矩阵，输出 `2x3` 子图热力图 PNG。

### 具体实现方式

- **`src/SteppingAction.cc`**
  - 在 `volName.contains("SiPM")` 命中时，直接读取 `copyNo` 作为 `SiPMBlockID`，调用 `EventAction::AddPhotonAtFaceBlock(...)`。

- **`include/EventAction.hh` + `src/EventAction.cc`**
  - 事件内新增 `fEventSiPMCounts`（`SiPMBlockID -> count`）容器。
  - 事件结束时从 `SiPMBlockID` 反解 `crystalId/face/j/k`，结合阵列 `fNy/fNz` 反解 `iy/iz`，写入事件明细 CSV。

- **`include/HistoManager.hh` + `src/HistoManager.cc`**
  - 新增 `FillPhotonFaceBlockEvent(...)` 接口与两个 ntuple：`PhotonFaceBlockEvent`、`PhotonFaceBlockTotal`。
  - 在写入事件明细时，同时在 `HistoManager` 内累加运行总计；`Save()` 阶段统一落盘 `PhotonFaceBlockTotal`。

### 阶段验证结果

- `cmake --build build --config Release` 编译通过。
- 运行宏后已创建：
  - `AnaEx01_nt_PhotonFaceBlockEvent*.csv`
  - `AnaEx01_nt_PhotonFaceBlockTotal*.csv`
- 脚本验证：
  - `python plot_sipm_6faces_heatmap.py --csv <...PhotonFaceBlockEvent_t0.csv> --out <.../sipm_6faces_heatmap_t0.png>`
  - 成功输出 6 面 `4x4` 热力图。

### 对前序开发的修正/延续说明

- 本次不是否定“六面几何放置”实现，而是将统计链路补齐到与六面几何一致（避免旧 Left/Right 统计口径与新几何不一致）。
- 旧 `PhotonLeft/PhotonRight` 与 `PhotonLRPerRod` 输出仍保留（兼容已有分析脚本）；新分析建议优先使用 `PhotonFaceBlockEvent/Total`。

---

## 增量更新：停用旧统计CSV，仅保留6面统计CSV（2026-02-21）

### 调整目标

- 按当前阶段需求，仅保留 6 面 SiPM 新增统计输出：
  - `PhotonFaceBlockEvent`
  - `PhotonFaceBlockTotal`
- 旧版统计 CSV（如 `PhotonLRPerRod`、`PhotonRight`、`PhotonLeft`、`PhotonDepthNtuple` 等）暂时停用。

### 实现方式

- 在 `src/HistoManager.cc` 的 `Book()` 中，将旧 photon 相关 ntuple 创建代码整体注释化（保留注释便于后续恢复）。
- 在 `src/EventAction.cc` 的 `EndOfEventAction()` 中，将旧 photon CSV 填充调用注释化，仅保留 6 面块级统计写入。
- 在 `src/HistoManager.cc` 中对旧 `FillPhoton*` 接口加“停用保护”（no-op/guard），避免误调用时写入旧 ntuple。

### 结果

- 6 面统计链路保持可用，输出不受影响。
- 旧统计 CSV 不再创建/写入（按当前需求停用）。

---

## 增量更新：BlockTotal 空数据根因修复（2026-02-22）

### 问题现象

- `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv` 仅部分线程文件有数据；
- `AnaEx01_nt_PhotonFaceBlockTotal*.csv` 只有表头，无数据行。

### 根因

- `RunAction::EndOfRunAction()` 仅在 `nofEvents == 0` 分支调用了 `fHistoManager->Save()`；
- 正常有事件的运行路径未调用 `Save()`，导致 `HistoManager::Save()` 中写入 `PhotonFaceBlockTotal` 的逻辑没有执行。

### 修复

- 文件：`src/RunAction.cc`
- 在 `EndOfRunAction()` 正常路径末尾补充：
  - `fHistoManager->Save();`

### 说明

- `t0~t7` 是线程编号，不是面编号；某些线程事件少或未命中 SiPM 时，`BlockEvent_tN` 为空是可能的。
- 面编号仍是 `Face=0..5`。

---

## 增量更新：移除 BlockTotal，改为 Python 合并线程（2026-02-22）

### 调整内容

- 按需求移除 `PhotonFaceBlockTotal` 输出，仅保留 `PhotonFaceBlockEvent`。
- 根目录脚本 `plot_sipm_6faces_heatmap.py` 新增线程文件合并能力，可直接读取结果目录中的
  `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv` 并自动聚合绘图。

### 使用方式

- 单文件绘图：
  - `python plot_sipm_6faces_heatmap.py --csv <path_to_csv> --out <png_path>`
- 线程文件自动合并绘图：
  - `python plot_sipm_6faces_heatmap.py --dir <results_folder> --out <png_path>`

---

## 增量更新：Histo9 风格批处理输出（2026-02-22）

### 目标

- 脚本行为对齐 `Histo9.py`：
  - 输入 `Results` 根目录后遍历全部配置子目录；
  - 在当前目录 `Output` 下自动创建时间戳输出目录；
  - 每个配置子目录输出“合并CSV + 聚合CSV + 6面热力图”。

### 实现

- 文件：`plot_sipm_6faces_heatmap.py`
- 默认用法：
  - `python plot_sipm_6faces_heatmap.py Results`
- 输出结构：
  - `Output/SiPM6_Output_<timestamp>/<config_name>/merged_event.csv`
  - `Output/SiPM6_Output_<timestamp>/<config_name>/merged_face_jk.csv`
  - `Output/SiPM6_Output_<timestamp>/<config_name>/sipm_6faces_heatmap.png`
- 行为说明：
  - `merged_event.csv`：合并 `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv` 原始事件行；
  - `merged_face_jk.csv`：按 `Face,j,k` 聚合后的计数结果；
  - 无可用 `PhotonFaceBlockEvent` 的配置会跳过并打印 `[SKIP]`。

---

## 增量更新：SiPM/分析输出路径说明（2026-02-22）

### Geant4 CSV 输出路径（Results）

- HistoManager 使用**相对路径**：`Results/<时间戳>_Nx..._Ny.../`，所有 ntuple CSV（含 `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv` 等）均写在该目录下。
- **实际落盘位置 = 程序启动时的当前工作目录 (CWD) + `Results/` + 上述子目录**。
- 常见情况：
  - 在**项目根目录**下执行 `build\Release\exampleB1.exe run4.mac` → CSV 在 **根目录的 `Results\...`**。
  - 在 **build\Release** 下双击运行或从该目录启动（如 Qt/IDE 工作目录设为 build）→ CSV 在 **`build\Results\...`**。
- 若在 Qt 可视化下运行且根目录与 build 下都未找到，请查看控制台打印的 “Output file is open in ...” 或 “Results\...” 路径，或在磁盘上搜索 `AnaEx01_nt_PhotonFaceBlockEvent` 或带日期的 `Results` 文件夹以定位。

### Python 分析脚本输出路径（Output）

- **Histo10_Cubic.py**（及同逻辑的 `plot_sipm_6faces_heatmap.py`）默认在**当前工作目录**下创建：
  - `Output/SiPM6_Output_<时间戳>/<配置名>/merged_event.csv`
  - `Output/SiPM6_Output_<时间戳>/<配置名>/merged_face_jk.csv`
  - `Output/SiPM6_Output_<时间戳>/<配置名>/sipm_6faces_heatmap.png`
- 使用 `--output <路径>` 可指定其它输出根目录。

---

## 增量更新：6面热力图坐标轴语义优化（2026-02-27）

### 调整目标

- 让 `Histo10_Cubic.py` 生成的 6 面热力图在每个子图上明确标注横纵轴物理含义，避免仅显示 `j/k` 时阅读歧义。

### 实现方式

- 文件：`Histo10_Cubic.py`
- 新增 `FACE_AXIS_LABELS` 映射，按 `face=0..5` 分别定义子图 `x/y` 轴标签。
- 保持数据矩阵填充逻辑不变（`matrix[j, k]`），仅增强坐标注释：
  - `imshow` 横轴仍对应 `k`，纵轴仍对应 `j`；
  - 轴标签中补充其在世界坐标中的方向（如 `k (+Z)`、`j (+Y)`）。

### 结果

- 每个面的子图都能直接看出“该面法向 + 本面内 j/k 对应的空间方向”，便于和 `DetectorConstruction` 的几何定义做人工比对。
- 未改变聚合口径与计数结果，不影响已有 `merged_event.csv` / `merged_face_jk.csv` 数据解释。

---

## 增量更新：3D 方位角交互工具性能优化（2026-03-06）

### 问题现象

- `cube_azimuth_interactive.py` 在拖动滑块时明显卡顿，尤其是连续拖动 `x/y/z` 或网格数滑块时。

### 根因分析

- `update()` 中每次交互都执行了大量“重建型”操作：
  - `ax3d.cla()` 后重画三面墙和立方体边线；
  - 三个热力图轴 `cla()` + `imshow()` 重建；
  - 三个 colorbar 通过 `cla()` + `fig.colorbar(...)` 重建；
  - 文本面板 `cla()` 后重建全部文本对象。
- 这些操作会在拖动过程中高频触发，造成 UI 线程负担过大。

### 实现方式

- 文件：`cube_azimuth_interactive.py`
- 采用“初始化一次 + 增量更新”策略：
  - 3D 静态元素（墙面、边线、坐标轴）只绘制一次；
  - 点源使用已存在的 `scatter` 对象，仅更新 `_offsets3d`；
  - 三个热力图只创建一次，交互时仅 `set_data()` + `set_extent()`；
  - colorbar 只创建一次，不再每帧重建；
  - 说明文本静态部分固定，仅更新位置字符串 `set_text()`；
  - 新增 `last_state` 缓存，避免相同滑块状态重复重绘。

### 结果与影响

- 交互响应明显提升，连续拖动时掉帧和迟滞显著减少。
- 功能保持不变：仍支持三面独立网格数、实时方位角计算、三幅热力图联动显示。
- 代码结构更清晰，便于后续扩展（例如增加极角、切换 colormap、导出当前矩阵等）。

---

## 增量更新：交互图中文字符缺失修复（2026-03-06）

### 问题现象

- `cube_azimuth_interactive.py` 的说明文本中，中文显示为方框或空白字符。

### 根因

- Matplotlib 默认字体在当前环境不一定包含中文字形，导致 CJK 字符无法渲染。

### 修复方式

- 文件：`cube_azimuth_interactive.py`
- 新增 `configure_plot_fonts()`，在程序启动时配置字体回退链：
  - `Microsoft YaHei`、`SimHei`、`Noto Sans CJK SC`、`Arial Unicode MS`、`DejaVu Sans`
- 同时设置 `axes.unicode_minus = False`，避免某些字体下负号显示异常。

### 结果

- 中文说明文本在 Windows 环境下显示稳定性明显提升；无对应字体时会自动回退到后备字体。

---

## 增量更新：方位角交互界面英文统一（2026-03-06）

### 调整内容

- 文件：`cube_azimuth_interactive.py`
- 将文本面板中的中文说明改为全英文：
  - “Azimuth definition (local wall coordinates) ...”
  - “Angle range: 0 deg ~ 360 deg”
  - “Current source position: (...)”
- 删除此前为中文渲染添加的 `configure_plot_fonts()` 及其调用，保持脚本更简洁。

### 结果

- 交互界面文本已统一为英文，避免中文字体兼容问题。
- 功能与性能优化逻辑保持不变。

---

## 增量更新：矩形精确立体角逻辑接入（2026-03-06）

### 调整目标

- 将交互脚本中的立体角计算从“未实现/近似思路”升级为“矩形四角点精确公式”。
- 对每个网格单元按边界坐标 `(x1,x2,y1,y2)` 计算立体角，避免中心点近似误差。

### 实现方式

- 文件：`cube_azimuth_interactive.py`
- 新增函数：
  - `rectangular_solid_angle(x1, x2, y1, y2, z)`  
    使用角点组合：
    `Omega = f(x2,y2)-f(x1,y2)-f(x2,y1)+f(x1,y1)`，
    `f(x,y)=atan2(x*y, z*sqrt(x^2+y^2+z^2))`
  - `wall_solid_angle(sx, sy, sz, wall, n)`  
    在每个墙面上把网格边界映射到局部坐标后，批量计算每个 cell 的 `Omega`。
- 在交互更新中实时计算三面 `Omega` 矩阵，并显示每面总立体角（sr）：
  - `x=0 sum(Omega_x0)`
  - `y=0 sum(Omega_y0)`
  - `z=0 sum(Omega_z0)`

### 结果

- 立体角逻辑与理论公式一致，`x_i/y_j` 按边界坐标组合实现，几何解释清晰。
- 保留原有方位角热力图显示，不改变原交互主视图。

---

## 增量更新：热力图改为立体角主视图，修正“中心点不对称”认知问题（2026-03-06）

### 问题说明

- 用户观察到点源接近中心时热力图看起来“不对称”。
- 核对后发现：
  - 立体角数值计算本身是正确的；
  - 但界面主图仍显示的是 `azimuth`，不是 `solid angle`；
  - 且截图中的点位置为 `(0.50, 0.49, 0.49)`，并非严格几何中心。

### 数值核对

- 对 `sx=sy=sz=0.5`，三面总立体角均为：
  - `2.0943951023931957 sr`
- 与立方体中心点按对称性应得到的单面立体角一致（`4π/6`）。

### 调整内容

- 文件：`cube_azimuth_interactive.py`
- 将三幅主热力图从 `azimuth` 改为 `solid angle`：
  - 初始化使用 `wall_solid_angle(...)`
  - 更新时显示 `Omega` 矩阵而不是 `azimuth` 矩阵
  - colorbar 单位改为 `sr`
  - 标题改为 `wall solid angle`
- 初始点位置改为严格中心：
  - `(0.50, 0.50, 0.50)`
- 文本说明改为矩形精确立体角公式说明，避免“图显示的是什么量”产生歧义。

### 结果

- 当点源位于中心时，三面热力图与总立体角表现符合对称直觉。
- 当前界面主视图与立体角逻辑保持一致，不再出现“后台算的是立体角，前台显示的是方位角”的错位。

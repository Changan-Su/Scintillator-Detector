# Log: 连续无空隙晶体条实现

**日期**: 2026-02-14  
**项目**: Scintillator-Detector-Continuous  
**版本**: v1.4

---

## 开发内容

将晶体阵列改为支持真正的无空隙连续排布（`crystalGap = 0 mm`），并在 gap 为 0 时自动禁用 Fillter，避免零厚度几何体。同时同步更新 `geometry.mac` 与批处理脚本默认参数，使快捷入口默认即为连续晶体配置。

---

## 实现方式

### 1) 放开 `crystalGap = 0` 参数校验

- 文件：`include/DetectorConstruction.hh`
- 修改 `SetCrystalGap` 的保护逻辑：
  - 旧逻辑：`v > 0` 才有效，否则回退到 `0.1`
  - 新逻辑：`v >= 0` 有效，仅负值回退到 `0.1`
- 结果：`geometry.mac` 中 `/detector/crystalGap 0 mm` 现在可真实生效。

### 2) gap=0 时禁用 Fillter 几何创建与放置

- 文件：`src/DetectorConstruction.cc`
- 关键改动：
  - `IfFillter` 改为 `Crystal_gap > 0.`
  - 仅在 `IfFillter` 为真时创建 `Fillter` 的 solid/logical volume
  - 放置条件改为 `ix != 0 && IfFillter && logicFillter != nullptr`
- 结果：无 gap 配置下不再出现 Fillter 体，不会产生零厚度 Fillter 的几何风险。

### 3) 同步宏与批处理默认参数

- 文件：`geometry.mac`
  - 保持 `/detector/crystalGap 0 mm`
  - 注释明确：`0` 表示连续无 gap 排布
  - 增加说明：gap=0 时 Fillter 自动禁用
- 文件：`run_batch.bat`
  - `GAP_START/GAP_END` 默认改为 `0.0`
  - `DEFAULT_GAP` 改为 `0.0`
  - 自动生成 `geometry.mac` 的注释改为连续无 gap 语义
- 文件：`run_batch.sh`
  - `GAP_START/GAP_END` 默认改为 `0.0`
  - `DEFAULT_GAP` 改为 `0.0`
  - 自动生成 `geometry.mac` 的注释改为连续无 gap 语义

---

## 验证记录

### 构建验证

- 发现原 `build` 目录缓存指向其他工程路径（`Scintillator-Detector-Single Rod`），会导致错误目标编译。
- 处理方式：新建并使用独立构建目录 `build_continuous`，并显式指定：
  - `Geant4_DIR=D:/Geant4/geant4-install/lib/cmake/Geant4`
  - `CMAKE_PREFIX_PATH=D:/Qt2/5.15.2/msvc2019_64/lib/cmake`
- 构建成功输出：`build_continuous/Release/exampleB1.exe`

### 运行验证

- 使用 `run4.mac` 进行 smoke test（100 events）。
- 结果：
  - 日志中不再出现 `Unit_fillter`
  - Overlap 检查均为 `OK`
  - `DEBUG:CRYSTAL_X` 显示 `16.5`（对应 `Nx=11, size=3mm, gap=0`，符合连续排布）
  - 输出 CSV 正常生成

---

## 对既有开发内容的修正说明

- 之前 `geometry.mac` 虽写 `0 mm`，但因 `SetCrystalGap()` 保护逻辑，实际会退回 `0.1 mm`，导致“无空隙”名义与实际不一致。
- 本次修正后，宏配置与真实几何行为一致，连续晶体条可复现。

---

## 增量更新（CSV 统一输出目录）

### 开发内容

- 将 CSV 输出从“依赖当前工作目录”改为统一输出到 `Results/时间+参数/`。
- 保持 CSV 文件名前缀不变，仍为 `AnaEx01_nt_*.csv`。
- 统一兼容两种入口：批处理宏运行与 Qt 可视化交互运行。

### 实现方式

- 文件：`src/HistoManager.cc`、`include/HistoManager.hh`
  - 在 `Book()` 阶段统一构建输出路径：`Results/<timestamp+geometry-params>/AnaEx01`。
  - 目录名包含 `Nx/Ny/Nz/Gap/Size/SizeY/FRY/FRZ/FPRY/FPRZ` 参数。
  - 增加毫秒时间戳与 `_dupN` 防冲突后缀，避免同秒重名。
  - 主/工线程共享同一输出路径，避免多线程写入分裂到多个目录。
- 文件：`include/DetectorConstruction.hh`
  - 补充 Fillter 参数 getter，供输出目录命名读取。
- 文件：`run_batch.bat`
  - 移除“根目录收集 CSV 后再搬运”的旧逻辑。
  - 批处理直接依赖程序写入 `Results/时间+参数/`，并打印本轮最新输出目录。

### 验证结果

- 从项目根目录执行 `build\\Release\\exampleB1.exe run1.mac` 后，CSV 已写入：
  - `Results/20260214_235057_148_Nx9_Ny1_Nz1_Gap0_Size3_SizeY3_FRY0p3_FRZ1_FPRY0p9_FPRZ0/`
- 根目录未再产生新的 `AnaEx01_nt_*.csv` 散落文件。

---

## 增量更新（Hosto9 连续晶体稳健峰检）

### 开发内容

- 新增 `Hosto9.py` 作为 `Histo8.py` 的连续晶体测量稳健版分析脚本。
- 针对“连续晶体下峰形异常/疑似假峰”问题，调整峰检策略以降低噪声峰误检。

### 实现方式

- 文件：`Hosto9.py`
  - 峰检测基于未加权计数直方图（count-based），不再用 `L+R` 加权曲线直接找峰。
  - 默认不强制目标峰数，不再通过逐步降低阈值去“凑峰”。
  - 增加事件门限参数：仅使用 `L+R >= min_total_photons` 的事件进行峰检（默认 `10`）。
  - 增加峰 prominence 绝对下限（默认 `5`），并保留相对阈值共同约束。
  - 输出图同时展示：
    - 归一化计数平滑曲线（用于找峰）
    - 归一化 `L+R` 加权平滑曲线（用于对比形状偏移）

### 验证结果

- 运行命令：
  - `python Hosto9.py Results --output Hosto9_Output_preview`
- 在当前两组可解析连续晶体数据上，峰数结果由原先易出现多峰/伪峰收敛为单峰：
  - `20260215_094445...`：`NumPeaks=1`，峰位约 `0.5117`
  - `20260215_102926...`：`NumPeaks=1`，峰位约 `0.8917`（低统计样本，仍建议增大事件数复核）
- 输出目录：
  - `Hosto9_Output_preview/`

### 对既有分析方式的修正说明

- 原 `Histo8.py` 在连续晶体工况中，强制目标峰数与降阈值策略可能将噪声起伏识别为真实峰。
- `Hosto9.py` 的定位是先保证峰检保守和稳定，再根据实验目标决定是否恢复多峰拟合策略。

---

## 增量更新（仅统一 Histo9.py 的 DOI 配置）

### 开发内容

- 用户仅保留并使用 `Histo9.py`，将 DOI 定义改为可配置模式，避免“改了公式但运行脚本不一致”。
- 将“事件门限/加权所用总光子数”与“DOI 横轴定义”解耦，避免 DOI 公式改动误伤筛选逻辑。

### 实现方式

- 文件：`Histo9.py`
  - 新增 DOI 模式参数 `--doi-mode`，支持：
    - `r_over_sum`：`R/(L+R)`
    - `r_over_l`：`R/L`
    - `diff_over_sum`：`(R-L)/(R+L)`
    - `log_r_over_l`：`log(R/L)`
  - 新增 `--doi-min` / `--doi-max`，并按 DOI 模式提供默认区间。
  - 固定使用 `L+R` 作为事件门限与加权基准，避免因 DOI 定义变化导致门限语义漂移。
  - 图标题与 x 轴标签自动显示当前 DOI 公式；输出 `summary.csv` 记录 `DOIMode/DOIRange`。
  - 输出目录前缀统一为 `Histo9_Output_*`（不再混用 Hosto9 命名）。

### 验证结果

- 命令：
  - `python Histo9.py Results --output Histo9_Output_smoke`
  - `python Histo9.py Results --output Histo9_Output_r_over_l --doi-mode r_over_l`
- 两次均成功完成；`summary.csv` 正常生成并带 DOI 模式信息。

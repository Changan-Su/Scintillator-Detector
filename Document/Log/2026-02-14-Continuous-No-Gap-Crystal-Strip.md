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

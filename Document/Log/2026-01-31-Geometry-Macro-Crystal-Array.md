# Log: 几何配置脚本控制晶体排布与填充

**日期**: 2026-01-31  
**项目**: Scintillator-Detector-Single Rod  
**版本**: v1.2

---

## 开发内容

通过**一个配置文件脚本**（Geant4 宏文件 `geometry.mac`）在 `/run/initialize` 之前设置晶体阵列（nx, ny, nz）和晶体间隙（Crystal_gap），以及可选单晶尺寸（crystalSize / crystalSizeY），无需修改 C++ 代码即可切换几何。同时支持在 Qt 可视化界面里交互修改参数并实时重建几何。

---

## 实现方式

### 1. 修改 `include/DetectorConstruction.hh`

- `fCrystal_gap` 类型由 `G4int` 改为 `G4double`（单位 mm），默认 0.1
- `fPar_ny`、`fPar_nz` 默认值改为 7，与原有 Construct() 行为一致
- `fcrystal_l`、`fcrystal_ly` 增加默认值 3.0（mm）
- 增加 setter：`SetCrystalGap(G4double)`、`SetCrystalSize(G4double)`、`SetCrystalSizeY(G4double)`，内部做 >0 保护
- `GetCrystal_gap()` 返回类型改为 `G4double`
- 构造函数/析构函数改为非 default，以便在构造函数中创建 DetectorMessenger、析构中释放
- 在 `B1` 命名空间内前向声明 `DetectorConstruction` 使用的 `DetectorMessenger` 类型

### 2. 修改 `src/DetectorConstruction.cc`

- 增加 `#include "DetectorMessenger.hh"`
- 实现构造函数：`fMessenger = new DetectorMessenger(this)`
- 实现析构函数：`delete fMessenger`
- 在 `Construct()` 中删除局部常量，改为使用成员变量：
  - `Crystal_gap = fCrystal_gap * mm`，`Crystal_nx/ny/nz = fPar_nx/ny/nz`，`crystal_l = fcrystal_l * mm`，`crystal_ly = fcrystal_ly * mm`
- 几何构建完成后，将长度写回成员时统一为 mm 数值（如 `fCrystal_gap = Crystal_gap / mm`）

### 3. 新增 `include/DetectorMessenger.hh` 与 `src/DetectorMessenger.cc`

- `DetectorMessenger` 继承 `G4UImessenger`，持有 `DetectorConstruction*`
- 命令目录：`/detector/`
- 命令：
  - `/detector/arrayNx`、`arrayNy`、`arrayNz`（整数）
  - `/detector/crystalGap`、`crystalSize`、`crystalSizeY`（带单位，默认 mm）
  - `/detector/update`（无参数）：调用 `G4RunManager::ReinitializeGeometry()` 重建几何，使参数修改立即生效
- `SetNewValue()` 中根据命令调用 `DetectorConstruction` 的对应 setter 或执行几何更新
- CMake 使用 `file(GLOB ... src/*.cc)`，新增 `DetectorMessenger.cc` 会自动加入编译，无需改 CMakeLists 源列表

### 4. 新增 `geometry.mac` 与 run 宏修改

- 新建 **geometry.mac**：示例设置 arrayNx/Ny/Nz、crystalGap、crystalSize（注释 crystalSizeY 和 update）
- **run2.mac**、**run1.mac** 在 `/run/initialize` 前增加 `/control/execute geometry.mac`
- **CMakeLists.txt**：将 `geometry.mac` 加入 `EXAMPLEB1_SCRIPTS` 及 POST_BUILD 复制列表；POST_BUILD 改为从 **PROJECT_SOURCE_DIR** 复制宏到 exe 目录（不再从 PROJECT_BINARY_DIR），使每次构建后 exe 目录使用项目源码中的最新 .mac

### 5. 修改 `run_vis.bat`（v1.2）

- 从**项目根目录**启动 exe：运行 exe 时不再执行 `cd build`，改为在项目根执行 `build\Release\exampleB1.exe`。
- 这样程序**当前工作目录 = 项目根**，`/control/execute geometry.mac` 和 `init_vis.mac` 会从**项目根**读取，即用户编辑的那份 `geometry.mac`。
- **效果**：改项目根目录下的 `geometry.mac` 后，用 **run_vis.bat** 启动，**无需重新构建**即可生效；仅需关闭程序、保存 .mac、再运行 run_vis.bat。

---

## 使用说明

### 方式一：通过 geometry.mac 预设参数（批处理 / 启动时配置）

**要让改 geometry.mac 不重新构建即生效**：必须用 **run_vis.bat** 从项目根目录启动（双击项目根下的 run_vis.bat）。此时程序工作目录为项目根，读的是项目根下的 `geometry.mac` 和 `init_vis.mac`。若从 `build\Release\` 直接双击 exe 或从 IDE 运行，工作目录多为 exe 所在目录，读的是那里的宏，需重新构建后才会更新。

1. 编辑**项目根目录**下的 **geometry.mac** 修改晶体阵列与填充，例如：
   ```geant4
   /detector/arrayNx 11
   /detector/arrayNy 7
   /detector/arrayNz 7
   /detector/crystalGap 0.1 mm
   /detector/crystalSize 3 mm
   ```
2. 运行前**必须**在 `/run/initialize` 之前执行几何配置，例如在 run 宏中：
   ```geant4
   /control/execute geometry.mac
   /run/initialize
   ```
3. 可维护多份几何宏（如 `geometry_5x5.mac`），在 run 宏里只 `execute` 其中一份再 `initialize`，实现"只改一个脚本"切换几何。

### 方式二：在 Qt 可视化界面里交互修改（实时调整）

1. 启动 Qt 可视化界面：推荐用 **run_vis.bat**（从项目根启动）；或双击 `build\Release\exampleB1.exe`（此时读 exe 目录下的宏）
2. 在命令行里输入参数修改命令，例如：
   ```geant4
   /detector/arrayNy 5
   /detector/crystalGap 0.2 mm
   ```
3. **必须**执行 `/detector/update` 命令，触发几何重建：
   ```geant4
   /detector/update
   ```
   终端会输出 "Geometry updated with current parameters."，此时可视化窗口中的几何会立即更新为新参数。

**注意**：
- `/detector/update` 在内部调用 `G4RunManager::ReinitializeGeometry()`，会销毁旧几何并用当前成员变量重新构建。
- 若在 `/run/initialize` **之前**改参数，**不需要** `/detector/update`（因为几何还没建）。
- 若在 `/run/initialize` **之后**改参数（如 Qt 里交互修改），**必须** `/detector/update` 才能生效。

### 小结：改 geometry.mac 后何时需重新构建

| 启动方式 | 工作目录 | 读到的 geometry.mac | 改 .mac 后是否需重新构建 |
|----------|----------|---------------------|---------------------------|
| **run_vis.bat**（项目根） | 项目根 | 项目根下的文件 | **否**，保存后重启程序即可 |
| 双击 build\Release\exampleB1.exe 或 IDE 运行 | exe 目录 | exe 目录下的副本 | **是**，需重新构建后才会更新 |

---

## 参考

- 计划文档：几何配置脚本控制晶体排布
- [[Geant4-and-Scintillator-Detector]] 项目规则

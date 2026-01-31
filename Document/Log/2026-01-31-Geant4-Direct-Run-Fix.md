# Log: Geant4 exampleB1 直接运行修复

**日期**: 2026-01-31  
**项目**: Scintillator-Detector-Single Rod  
**版本**: v1.0

---

## 开发内容

修复 Geant4 Qt 可视化程序 `exampleB1.exe` 无法直接启动（双击或命令行运行）的问题，实现无需每次通过 `run_vis.bat` 即可直接运行。

---

## 实现方式

### 1. 新增 `CopyDlls.cmake` 脚本

- 在 POST_BUILD 阶段自动将 **Geant4 DLL**（G4*.dll、msvcp、vcruntime、concrt 等）复制到可执行文件目录
- 将 **Qt DLL**（Qt5*.dll 及 platforms 插件）复制到可执行文件目录
- 使用方式：`cmake -DDEST_DIR=<path> -DG4_BIN=<path> -DQT_BIN=<path> -P CopyDlls.cmake`

### 2. 修改 `CMakeLists.txt`

- 增加宏文件复制到 exe 目录的 POST_BUILD 步骤（init_vis.mac, vis.mac, run1.mac, run2.mac）
- 增加调用 `CopyDlls.cmake` 的 POST_BUILD 步骤，在每次构建后自动部署 DLL
- 通过 `Geant4_DIR` 和 `Qt5Core_DIR` 自动解析 Geant4、Qt 的 bin 路径

### 3. 新增 `run_vis.bat`（备用启动方式）

- 显式设置 Geant4、Qt 的 PATH 后启动 exampleB1
- 支持两种输出路径：
  - `build\Release\`（标准 CMake 构建）
  - `D:\Geant4\geant4-install\share\Geant4\examples\basic\B1\Release\`（旧 vcxproj 输出）

---

## 使用说明

### 直接运行（推荐）

1. 在项目根目录创建 build 并构建：
   ```bat
   mkdir build
   cd build
   cmake ..
   cmake --build . --config Release
   ```
2. 构建完成后，可直接双击 `build\Release\exampleB1.exe` 或在任意终端运行，无需 `run_vis.bat`

### 备用：通过 run_vis.bat 启动

当 DLL 未正确部署时，可使用 `run_vis.bat` 通过 PATH 加载依赖。

---

## 同步项目

相同修复已应用于 **Scintillator-Detector-master** 项目。

---

## 参考

- [[Geant4-and-Scintillator-Detector]] 项目规则
- Geant4 11.3.1 / Qt 5.15.2 msvc2019_64

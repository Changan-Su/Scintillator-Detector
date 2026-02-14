# Log: G4Tree.dll 找不到 — CopyDlls 路径修复

**日期**: 2026-02-14  
**项目**: Scintillator-Detector-Continuous  
**版本**: v1.5

---

## 开发内容

修复运行 `exampleB1.exe` 时提示「找不到 G4Tree.dll」的问题。根因是 `CMakeLists.txt` 中由 `Geant4_DIR` 推导 Geant4 `bin` 目录的层级错误，导致 POST_BUILD 的 `CopyDlls.cmake` 从未向可执行目录复制任何 Geant4 DLL（包括 G4Tree.dll）。

---

## 实现方式

### 1) 修正 GEANT4_BIN 推导路径

- **文件**: `CMakeLists.txt`
- **修改**:  
  `get_filename_component(GEANT4_BIN "${Geant4_DIR}/../../bin" ABSOLUTE)`  
  →  
  `get_filename_component(GEANT4_BIN "${Geant4_DIR}/../../../bin" ABSOLUTE)`
- **原因**:  
  - `Geant4_DIR` 实际为 `.../lib/cmake/Geant4`（例如 `D:/Geant4/geant4-install/lib/cmake/Geant4`）  
  - 原式 `../../bin` 得到的是 `.../lib/bin`（不存在）  
  - 正确应为从 `lib/cmake/Geant4` 上溯到安装根再进 `bin`，即 `../../../bin`

### 2) 验证

- 重新配置并构建：
  ```powershell
  cmake -S . -B build -DGeant4_DIR=D:/Geant4/geant4-install/lib/cmake/Geant4 -DCMAKE_PREFIX_PATH=D:/Qt2/5.15.2/msvc2019_64/lib/cmake
  cmake --build build --config Release
  ```
- 构建输出中应出现：`Copied 43 Geant4 DLLs to .../build/Release`
- `build/Release/G4Tree.dll` 及其余 G4*.dll 应存在，可直接运行 `exampleB1.exe`（双击或命令行）无需再设 PATH

---

## 对既有开发内容的修正说明

- 2026-01-31 的「Geant4 exampleB1 直接运行修复」中引入了 `CopyDlls.cmake` 与 `GEANT4_BIN` 推导；当时使用的 `Geant4_DIR/../../bin` 在**本机**安装结构（`Geant4_DIR = .../lib/cmake/Geant4`）下会指向错误目录，导致 CopyDlls 未复制任何 Geant4 DLL。本次修正后，与当前 Geant4 安装布局一致，DLL 部署正常。

---

## 参考资料

- 问题与排查过程见：`Document/Bugs/2026-02-14-G4Tree-DLL-NotFound.md`

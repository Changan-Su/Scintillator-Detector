# Bug Review: 运行时找不到 G4Tree.dll

**日期**: 2026-02-14  
**项目**: Scintillator-Detector-Continuous  
**Log 版本**: 2026-02-14-G4Tree-DLL-CopyDlls-Path-Fix

---

## 问题描述

- **现象**: 运行 `build/Release/exampleB1.exe`（双击或命令行）时，系统提示找不到 `G4Tree.dll`，程序无法启动。
- **环境**: Windows 10/11，Visual Studio 2022，CMake 生成 VS 工程，Geant4 安装于 `D:\Geant4\geant4-install`（例如 `Geant4_DIR = .../lib/cmake/Geant4`）。

---

## 原因分析

1. **DLL 未拷贝到 exe 目录**  
   项目通过 POST_BUILD 调用 `CopyDlls.cmake` 将 Geant4、Qt 的 DLL 复制到 `exampleB1.exe` 所在目录。检查发现 `build/Release/` 下没有任何 `G4*.dll`，说明 CopyDlls 使用的 Geant4 bin 路径错误，未找到源 DLL。

2. **GEANT4_BIN 推导错误**  
   - `CMakeLists.txt` 中原式为：  
     `get_filename_component(GEANT4_BIN "${Geant4_DIR}/../../bin" ABSOLUTE)`  
   - 本机 `Geant4_DIR` 为 `D:/Geant4/geant4-install/lib/cmake/Geant4`，则：  
     - `Geant4_DIR/../` = `.../lib/cmake`  
     - `Geant4_DIR/../../` = `.../lib`  
     - `Geant4_DIR/../../bin` = `.../lib/bin`（该目录不存在）  
   - 实际 Geant4 DLL 在 `.../geant4-install/bin`，即需要从 `lib/cmake/Geant4` 再向上一级到安装根，再进 `bin`，故应为 `Geant4_DIR/../../../bin`。

3. **PATH 未包含 Geant4**  
   若未将 `Geant4\bin` 加入系统或用户 PATH，且 exe 目录也没有拷贝 DLL，则运行时必然找不到 G4Tree.dll。本问题以「修正拷贝路径、保证 exe 旁有 DLL」为主，PATH 为可选备用。

---

## 解决方案

### 永久修复（已采用）

- **文件**: `CMakeLists.txt`
- **修改**: 将  
  `get_filename_component(GEANT4_BIN "${Geant4_DIR}/../../bin" ABSOLUTE)`  
  改为  
  `get_filename_component(GEANT4_BIN "${Geant4_DIR}/../../../bin" ABSOLUTE)`
- **操作**: 重新配置并构建后，POST_BUILD 会向 `build/Release` 复制 43 个 Geant4 DLL（含 G4Tree.dll），可直接运行 `exampleB1.exe`。

### 临时修复（不改代码）

若暂不修改 CMake，可先设置 PATH 再运行 exe：
```powershell
$env:PATH="D:\Geant4\geant4-install\bin;D:\Qt2\5.15.2\msvc2019_64\bin;$env:PATH"
.\build\Release\exampleB1.exe
```

---

## 验证

- 修改后执行：  
  `cmake -S . -B build -DGeant4_DIR=... -DCMAKE_PREFIX_PATH=...`  
  `cmake --build build --config Release`
- 检查：  
  `build/Release/G4Tree.dll` 存在，且可双击或命令行直接运行 `exampleB1.exe` 无报错。

---

## 参考资料

- 项目 Log：`Document/Log/2026-02-14-G4Tree-DLL-CopyDlls-Path-Fix.md`
- 直接运行方案（CopyDlls 引入）：`Document/Log/2026-01-31-Geant4-Direct-Run-Fix.md`

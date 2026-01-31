# Log: 参数扫描批量运行脚本

**日期**: 2026-01-31  
**项目**: Scintillator-Detector-Single Rod  
**版本**: v1.0

---

## 开发内容

创建参数扫描批量运行脚本（Windows `run_batch.bat` 和 Linux `run_batch.sh`），支持自定义几何参数循环（arrayNx/Ny/Nz、crystalGap、crystalSize、crystalSizeY），自动运行模拟并将结果按循环次数+参数命名存储到 `Results/` 目录，每个结果文件夹包含当次运行的 `geometry.mac` 副本。

---

## 实现方式

### 1. 脚本文件

- **run_batch.bat**（Windows 批处理脚本）
- **run_batch.sh**（Linux/Mac bash 脚本）
- **Results/.gitkeep**（确保 Results 目录被 git 追踪）

### 2. 配置区（脚本开头可编辑）

#### 参数循环配置（每个参数可独立开启/关闭）

| 配置项 | 说明 | 示例 |
|--------|------|------|
| `LOOP_ARRAY_NX` | 是否循环 arrayNx | `true` / `false` |
| `NX_START`, `NX_END`, `NX_STEP` | arrayNx 起始值、终止值、步长 | `11`, `15`, `2` |
| `LOOP_ARRAY_NY` | 是否循环 arrayNy | `true` / `false` |
| `NY_START`, `NY_END`, `NY_STEP` | arrayNy 起始值、终止值、步长 | `3`, `7`, `2` |
| `LOOP_ARRAY_NZ` | 是否循环 arrayNz | `true` / `false` |
| `NZ_START`, `NZ_END`, `NZ_STEP` | arrayNz 起始值、终止值、步长 | `3`, `7`, `2` |
| `LOOP_CRYSTAL_GAP` | 是否循环 crystalGap | `true` / `false` |
| `GAP_START`, `GAP_END`, `GAP_STEP` | crystalGap 起始值、终止值、步长（mm） | `0.1`, `0.5`, `0.1` |
| `LOOP_CRYSTAL_SIZE` | 是否循环 crystalSize | `true` / `false` |
| `SIZE_START`, `SIZE_END`, `SIZE_STEP` | crystalSize 起始值、终止值、步长（mm） | `2`, `4`, `1` |
| `LOOP_CRYSTAL_SIZE_Y` | 是否循环 crystalSizeY | `true` / `false` |
| `SIZE_Y_START`, `SIZE_Y_END`, `SIZE_Y_STEP` | crystalSizeY 起始值、终止值、步长（mm） | `2`, `4`, `1` |

#### 默认值配置（不循环时使用的固定值）

| 配置项 | 说明 | 默认值 |
|--------|------|--------|
| `DEFAULT_NX` | arrayNx 默认值 | `11` |
| `DEFAULT_NY` | arrayNy 默认值 | `7` |
| `DEFAULT_NZ` | arrayNz 默认值 | `7` |
| `DEFAULT_GAP` | crystalGap 默认值（mm） | `0.1` |
| `DEFAULT_SIZE` | crystalSize 默认值（mm） | `3` |
| `DEFAULT_SIZE_Y` | crystalSizeY 默认值（mm） | `3` |

#### 运行配置

| 配置项 | 说明 | 默认值 |
|--------|------|--------|
| `RUN_MACRO` | 使用的宏文件 | `run4.mac` |
| `EXE_PATH` | 可执行文件路径 | `build/Release/exampleB1.exe` (Windows) <br> `build/Release/exampleB1` (Linux) |

#### 文件夹命名配置（控制哪些参数显示在文件夹名中）

| 配置项 | 说明 | 默认值 |
|--------|------|--------|
| `NAME_INCLUDE_NX` | 文件夹名是否包含 Nx | `false` |
| `NAME_INCLUDE_NY` | 文件夹名是否包含 Ny | `false` |
| `NAME_INCLUDE_NZ` | 文件夹名是否包含 Nz | `false` |
| `NAME_INCLUDE_GAP` | 文件夹名是否包含 Gap | `true` |
| `NAME_INCLUDE_SIZE` | 文件夹名是否包含 Size | `false` |
| `NAME_INCLUDE_SIZE_Y` | 文件夹名是否包含 SizeY | `false` |

### 3. 脚本执行流程

1. 读取配置区，生成参数列表（仅循环开启的参数，其余使用默认值）
2. 创建 `Results/` 目录
3. 嵌套循环遍历所有参数组合：
   - 生成当前参数的 `geometry.mac`（覆盖项目根目录的 geometry.mac）
   - 创建结果文件夹：`Results/NNN_Param1Value1_Param2Value2_...`（NNN 为 3 位数循环序号，参数名根据 `NAME_INCLUDE_*` 配置）
   - 复制 `geometry.mac` 到结果文件夹
   - 运行 `exampleB1.exe RUN_MACRO`
   - 移动生成的 `AnaEx01_nt_*.csv` 到结果文件夹
4. 输出完成信息及总运行次数

### 4. 输出文件结构

```
Scintillator-Detector-Single Rod/
  Results/
    001_Gap0.1/
      AnaEx01_nt_Ntuple1.csv
      AnaEx01_nt_Ntuple2.csv
      AnaEx01_nt_PhotonLeft.csv
      AnaEx01_nt_PhotonRight.csv
      ... (其他 CSV 文件)
      geometry.mac              # 当次运行的参数配置副本
    002_Gap0.2/
      ...
    003_Gap0.3/
      ...
```

---

## 使用说明

### Windows

1. 编辑 **run_batch.bat** 配置区，设置要循环的参数及范围
2. 保存并双击运行 `run_batch.bat`
3. 等待批量运行完成（终端会显示进度）
4. 结果保存在 `Results/` 目录

### Linux/Mac

1. 编辑 **run_batch.sh** 配置区，设置要循环的参数及范围
2. 如果需要，在脚本中取消注释并调整 Geant4 环境变量加载路径：
   ```bash
   source /path/to/geant4-install/bin/geant4.sh
   ```
3. 添加执行权限并运行：
   ```bash
   chmod +x run_batch.sh
   ./run_batch.sh
   ```
4. 结果保存在 `Results/` 目录

### 配置示例一：扫描 crystalGap 从 0.1 到 0.5，步长 0.1

**Windows (run_batch.bat)**:
```bat
set LOOP_CRYSTAL_GAP=true
set GAP_START=0.1
set GAP_END=0.5
set GAP_STEP=0.1
set NAME_INCLUDE_GAP=true

REM 其他参数设为 false
set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NY=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false
```

**结果文件夹**（5 次运行）：
- `001_Gap0.1`
- `002_Gap0.2`
- `003_Gap0.3`
- `004_Gap0.4`
- `005_Gap0.5`

### 配置示例二：同时扫描 arrayNy 和 crystalGap

**Windows (run_batch.bat)**:
```bat
set LOOP_ARRAY_NY=true
set NY_START=3
set NY_END=7
set NY_STEP=2

set LOOP_CRYSTAL_GAP=true
set GAP_START=0.1
set GAP_END=0.3
set GAP_STEP=0.1

set NAME_INCLUDE_NY=true
set NAME_INCLUDE_GAP=true
```

**结果文件夹**（3×3 = 9 次运行）：
- `001_Ny3_Gap0.1`, `002_Ny3_Gap0.2`, `003_Ny3_Gap0.3`
- `004_Ny5_Gap0.1`, `005_Ny5_Gap0.2`, `006_Ny5_Gap0.3`
- `007_Ny7_Gap0.1`, `008_Ny7_Gap0.2`, `009_Ny7_Gap0.3`

### 配置示例三：只改一次 arrayNx，不循环

**Windows (run_batch.bat)**:
```bat
REM 所有循环都设为 false
set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NY=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_GAP=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false

REM 修改默认值
set DEFAULT_NX=15
set DEFAULT_NY=7
set DEFAULT_NZ=7

set NAME_INCLUDE_NX=true
set NAME_INCLUDE_NY=true
set NAME_INCLUDE_NZ=true
```

**结果文件夹**（1 次运行）：
- `001_Nx15_Ny7_Nz7`

---

## 注意事项

1. **工作目录**：脚本从项目根目录运行，exe 的工作目录设为项目根，这样 `geometry.mac` 和 `run*.mac` 都从项目根读取
2. **环境变量**：
   - Windows 版本会自动调用 `D:\Geant4\geant4-install\bin\geant4.bat` 设置环境（路径硬编码在脚本中，若安装位置不同需手动修改）
   - Linux 版本需要在脚本中取消注释并配置 Geant4 环境变量加载路径
3. **CSV 文件移动**：脚本会移动所有 `AnaEx01_nt_*.csv` 文件到对应结果文件夹；如果之前运行留下了旧 CSV，会一并移动（建议每次批量运行前清理项目根目录的 CSV）
4. **geometry.mac 覆盖**：脚本会在每次循环时覆盖项目根目录的 `geometry.mac`；批量运行完成后，`geometry.mac` 保持最后一次循环的参数
5. **步长精度**：脚本支持一位小数的步长（如 0.1、0.2），更高精度需修改脚本中的浮点数处理逻辑

---

## 参考

- [[2026-01-31-Geometry-Macro-Crystal-Array]] 几何配置脚本控制晶体排布与填充
- [[Geant4-and-Scintillator-Detector]] 项目规则

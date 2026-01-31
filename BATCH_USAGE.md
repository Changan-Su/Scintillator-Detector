# 参数扫描批量运行脚本使用示例

本文档演示如何使用 `run_batch.bat` (Windows) 或 `run_batch.sh` (Linux) 进行参数扫描。

## 快速开始

### Windows

1. 打开 `run_batch.bat`
2. 编辑配置区（脚本开头的 `CONFIGURATION SECTION`）
3. 保存并双击运行 `run_batch.bat`
4. 结果自动保存到 `Results/` 目录

### Linux/Mac

1. 打开 `run_batch.sh`
2. 编辑配置区（脚本开头的 `CONFIGURATION SECTION`）
3. 如果需要，配置 Geant4 环境变量路径
4. 运行：
   ```bash
   chmod +x run_batch.sh
   ./run_batch.sh
   ```

## 示例 1：扫描晶体间隙 (crystalGap)

扫描 0.1 mm 到 0.5 mm，步长 0.1 mm（共 5 次运行）

```bat
REM Windows 配置
set LOOP_CRYSTAL_GAP=true
set GAP_START=0.1
set GAP_END=0.5
set GAP_STEP=0.1
set NAME_INCLUDE_GAP=true

set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NY=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false
```

**结果**：
```
Results/
  001_Gap0.1/
  002_Gap0.2/
  003_Gap0.3/
  004_Gap0.4/
  005_Gap0.5/
```

## 示例 2：扫描阵列大小 (arrayNy)

扫描 y 方向晶体数 3, 5, 7（共 3 次运行）

```bat
set LOOP_ARRAY_NY=true
set NY_START=3
set NY_END=7
set NY_STEP=2
set NAME_INCLUDE_NY=true

set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_GAP=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false
```

**结果**：
```
Results/
  001_Ny3/
  002_Ny5/
  003_Ny7/
```

## 示例 3：二维参数扫描

同时扫描 arrayNy (3, 5, 7) 和 crystalGap (0.1, 0.2, 0.3)（共 3×3 = 9 次运行）

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

REM 其他参数关闭
set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false
```

**结果**：
```
Results/
  001_Ny3_Gap0.1/
  002_Ny3_Gap0.2/
  003_Ny3_Gap0.3/
  004_Ny5_Gap0.1/
  005_Ny5_Gap0.2/
  006_Ny5_Gap0.3/
  007_Ny7_Gap0.1/
  008_Ny7_Gap0.2/
  009_Ny7_Gap0.3/
```

## 示例 4：单次运行（指定特定参数）

不循环任何参数，只用指定的默认值运行一次

```bat
set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NY=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_GAP=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false

REM 设置要使用的参数
set DEFAULT_NX=15
set DEFAULT_NY=9
set DEFAULT_NZ=9
set DEFAULT_GAP=0.15
set DEFAULT_SIZE=4

set NAME_INCLUDE_NX=true
set NAME_INCLUDE_NY=true
set NAME_INCLUDE_NZ=true
set NAME_INCLUDE_GAP=true
set NAME_INCLUDE_SIZE=true
```

**结果**：
```
Results/
  001_Nx15_Ny9_Nz9_Gap0.15_Size4/
```

## 每个结果文件夹包含

- **AnaEx01_nt_*.csv**：Geant4 模拟生成的所有 CSV 数据文件
- **geometry.mac**：当次运行使用的几何参数配置副本

## 注意事项

1. **运行前清理**：建议删除项目根目录下的旧 CSV 文件，避免误移动到结果文件夹
2. **RUN_MACRO 配置**：默认使用 `run4.mac`（100000 events），可修改为 `run3.mac`（500 events）用于快速测试
3. **Geant4 环境**：
   - Windows 脚本会自动加载 `D:\Geant4\geant4-install\bin\geant4.bat`
   - Linux 脚本需手动配置 Geant4 环境变量路径（在脚本中取消注释相应行）
4. **步长精度**：当前脚本支持一位小数步长（如 0.1），需要更高精度请修改脚本

## 查看结果

每个文件夹中的 CSV 可用 Python 分析脚本处理，例如：

```bash
cd Results/001_Gap0.1
python ../../Histo7.py
```

或批量处理所有结果：

```bash
for dir in Results/*/; do
    cd "$dir"
    python ../../Histo7.py
    cd ../..
done
```

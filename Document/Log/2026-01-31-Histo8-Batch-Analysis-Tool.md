# 2026-01-31 Histo8 批量数据分析工具开发

## 版本信息
- **日期**: 2026-01-31
- **最终版本**: v1.1
- **相关文件**: 
  - `Histo8.py` (新增并修复)
  - `Histo7.py` (修复)

## 版本历史
- **v1.0**: 初始版本，基本功能实现
- **v1.1**: 修复 CSV 读取问题，改进 Summary 图为叠加显示

---

## 开发背景

在参数扫描批量运行脚本（`run_batch.bat`）生成大量配置数据后，需要一个工具来批量分析所有配置下每个 rod 的闪烁位置分布（DOI - Depth of Interaction）。原有的 `Histo7.py` 只能单个 rod 手动分析，效率低下。

## 主要开发内容

### 1. Histo8.py - 批量分析工具（新增）

#### 核心功能

**输入处理**：
- 自动扫描 `Results/` 文件夹下所有配置子文件夹
- 解析每个文件夹中的 `geometry.mac` 提取晶体配置参数
  - `arrayNx`, `arrayNy`, `arrayNz` - 晶体阵列尺寸
  - `crystalGap` - 晶体间隙
  - `crystalSize`, `crystalSizeY` - 晶体尺寸
- 自动检测并合并多线程 CSV 文件
  - 如果只有 `AnaEx01_nt_PhotonLRPerRod_t*.csv`，自动合并为单个文件
  - 使用 pandas concat 合并所有线程数据

**数据分析**：
- 基于 Histo7 的分析核心，对每个 rod 执行：
  1. 计算 DOI 不对称度：`asym = R / (L + R)`
  2. 生成加权直方图（权重为 L+R）
  3. 高斯平滑处理
  4. 自适应峰检测（多个 prominence 候选值）
  5. 多峰高斯拟合
  6. 计算统计指标：
     - 峰数量 (NumPeaks)
     - 平均半峰宽 (AvgFWHM)
     - 最小峰谷比 (MinPVR) 及分贝值

**输出生成**：
1. **单个 rod 直方图** (`rod_iy*_iz*.png`)
   - 显示原始直方图、平滑曲线、拟合曲线
   - 标注峰位置
   - 显示配置信息和统计指标
   - 150 DPI PNG 格式

2. **配置汇总图** (`summary.png`)
   - 网格化展示该配置所有 rod 的直方图
   - 使用 GridSpec 布局，按 (iy, iz) 排列
   - 每个子图显示峰数和 FWHM

3. **全局汇总 CSV** (`summary.csv`)
   - 包含所有配置、所有 rod 的统计数据
   - 字段：Config, iy, iz, Nx, Ny, Nz, Gap_mm, Size_mm, SizeY_mm, NumEvents, NumPeaks, AvgFWHM, MinPVR, MinPVR_dB

#### 实现方式

**主要函数**：

1. `parse_geometry_mac(mac_path)`: 解析 geometry.mac 文件
   - 使用正则表达式匹配 `/detector/参数名 值` 格式
   - 返回配置字典

2. `merge_thread_csvs(folder, pattern)`: 自动合并多线程 CSV
   - 检查是否已有合并文件
   - 查找所有 `_t*.csv` 文件
   - 使用 pandas concat 合并
   - 保存为单个 CSV

3. `read_photon_lr(csv_path)`: 鲁棒 CSV 读取
   - 优先尝试带表头读取
   - 自动大小写匹配列名
   - 回退到无表头模式并强制列名

4. `analyze_one_rod(df, iy, iz)`: 单个 rod 分析
   - 完整的直方图分析流程
   - 返回统计结果字典
   - 包含成功/失败状态和错误信息

5. `plot_single_rod(result, config, config_name, output_path)`: 绘制单个 rod 图
   - 使用 matplotlib 生成标注完整的直方图
   - 标题包含配置信息和统计指标

6. `plot_summary_grid(results, config, config_name, output_path)`: 绘制汇总网格图
   - 动态计算网格大小（ny × nz）
   - 每个子图显示简化的直方图信息

7. `process_config_folder(folder, output_folder)`: 处理单个配置
   - 解析配置、读取数据、分析所有 rod
   - 生成所有输出文件
   - 返回统计结果列表

**命令行接口**：
```bash
python Histo8.py D:\path\to\Results
python Histo8.py --results D:\path\to\Results --output ./output
```

**输出目录结构**：
```
Histo8_Output_YYYYMMDD_HHMMSS/
├── 001_Ny1/
│   ├── rod_iy0_iz0.png
│   ├── rod_iy0_iz1.png
│   └── summary.png
├── 002_Ny2/
│   └── ...
└── summary.csv
```

#### 参数配置

脚本顶部参数区（可按需调整）：
```python
n_bins = 300                 # DOI 直方图的 bin 数
sigma_smooth = 1.0           # 高斯平滑 sigma
distance_ratio = 0.035       # 峰间最小距离比例
prominence_ratio = 0.01      # 初始 prominence 比例
target_n = 11                # 期望峰数
sigma_init_guess = 0.010     # 高斯初始 sigma
mu_window = 0.030            # μ 搜索窗口
sigma_bounds = (0.003, 0.050)# sigma 上下界
```

### 2. Histo7.py 问题修复

修复了以下 bug：

#### Bug 1: 文档注释脚本名称错误
**问题**：文档注释中的用法示例仍引用 `Histo6.py`
```python
# 错误
python Histo6.py --iy 5 --iz 7
```
**修复**：更新为正确的脚本名称
```python
# 正确
python Histo7.py --iy 5 --iz 7
```

#### Bug 2: 第二次拟合缺少异常处理
**问题**：第二次 `curve_fit` 调用（放宽 sigma 上界重试）没有异常处理，如果再次失败会导致程序崩溃
```python
# 原代码
except Exception:
    ub2 = ub[:]
    for i in range(2, len(ub2), 3):
        ub2[i] = max(0.08, ub2[i])
    popt, _ = curve_fit(...)  # 如果失败会崩溃
    fitted_y = multi_gaussian(hist_x, *popt)
```

**修复**：添加嵌套 try-except，失败时输出友好提示
```python
# 修复后
except Exception as e1:
    try:
        ub2 = ub[:]
        for i in range(2, len(ub2), 3):
            ub2[i] = max(0.08, ub2[i])
        popt, _ = curve_fit(...)
        fitted_y = multi_gaussian(hist_x, *popt)
    except Exception as e2:
        print(f"Warning: 拟合失败 - {str(e2)}")
        print("将继续显示未拟合的直方图")
```

#### Bug 3: 注释缩进错误
**问题**：第 169-170 行注释缩进混乱
```python
    plt.show()

    # 打印拟合结果（可选）
        # 打印拟合结果 + 平均 FWHM
    if popt is not None:
```

**修复**：统一缩进格式
```python
    plt.show()

    # 打印拟合结果 + 平均 FWHM
    if popt is not None:
```

### 3. 文档创建

#### HISTO8_USAGE.md
创建了详细的使用指南文档，包含：
- 功能概述
- 使用方法和命令示例
- 输出结构说明
- 输出文件详细解释
- 数据要求
- 故障排除指南
- 与 Histo7.py 的对比
- 依赖库说明
- 性能提示

## 技术细节

### 关键技术点

1. **鲁棒的 CSV 读取**
   - 兼容带/不带表头的 CSV
   - 自动大小写匹配列名
   - 处理分隔符不一致问题

2. **自适应峰检测**
   - 多个 prominence 候选值递减尝试
   - 根据 target_n 筛选最显著的峰
   - 使用 scipy.signal.find_peaks 的完整参数

3. **双重拟合策略**
   - 首次拟合使用严格的 sigma 上界
   - 失败后放宽 sigma 上界重试
   - 两次都失败则跳过拟合但继续流程

4. **批量处理流程**
   - 文件夹遍历 → 配置解析 → 数据合并 → rod 分析 → 图表生成 → CSV 汇总
   - 单个 rod 失败不影响其他 rod
   - 所有输出带时间戳避免覆盖

### 依赖库
- `numpy`: 数值计算
- `pandas`: 数据处理和 CSV 操作
- `matplotlib`: 图表绘制
- `scipy`: 信号处理（平滑、峰检测）和拟合

## 使用示例

### 批量分析
```bash
# 分析 Results 文件夹中所有配置
python Histo8.py D:\Geant4\Projects\Scintillator-Detector-Single Rod\Results
```

### 单个 rod 分析（使用 Histo7）
```bash
# 分析某个配置的特定 rod
cd Results/001_Ny1
python ../../Histo7.py --iy 0 --iz 5 --csv AnaEx01_nt_PhotonLRPerRod.csv
```

## 验证测试

执行了以下验证：
1. ✅ Python 语法编译检查（`python -m py_compile`）
2. ✅ Histo8.py 无语法错误
3. ✅ Histo7.py 无语法错误
4. ✅ 修复了转义字符警告

## 性能考虑

- 使用 150 DPI 图片平衡质量与文件大小
- 批量处理大量配置可能耗时较长（取决于 rod 数量和事件数）
- 建议先用少量配置测试参数设置

## v1.1 更新内容（实际运行测试后的修复）

### Bug 修复 1: CSV 读取问题

**问题描述**：
初次运行时遇到 "Could not determine delimiter" 错误，原因：
1. Geant4 输出的 CSV 文件带有 `#` 开头的注释行
2. 合并后的 CSV 文件只有注释没有数据（空文件）
3. pandas 的 `sep=None, engine="python"` 无法从空文件推断分隔符
4. 线程文件使用逗号分隔符，但读取逻辑不够鲁棒

**修复方案**：
1. **优化 `read_photon_lr()` 函数**：
   - 优先使用逗号分隔符读取（Geant4 默认格式）
   - 正确处理 `comment="#"` 参数跳过注释行
   - 多重回退策略：逗号分隔 → 带表头 → 自动分隔符
   - 添加详细的错误提示

2. **改进 `merge_thread_csvs()` 函数**：
   - 检测空的合并文件并自动重新生成
   - 使用明确的逗号分隔符读取线程文件
   - 合并后保存为标准 CSV（带列名，无注释）

**修复后的代码逻辑**：
```python
# read_photon_lr - 优先使用逗号分隔符
df0 = pd.read_csv(csv_path, comment="#", header=None, sep=",")
df0 = df0.iloc[:, :5].copy()
df0.columns = ["EventID","iz","iy","Left","Right"]

# merge_thread_csvs - 检测空文件
with open(merged_path, 'r') as f:
    has_data = any(line.strip() and not line.startswith('#') 
                  for line in f.readlines())
if not has_data:
    # 重新合并
```

**测试结果**：
- ✅ 成功读取所有 3 个配置文件夹
- ✅ 自动合并 8 个线程文件（每个配置）
- ✅ 总计分析 18 个 rod，无错误

### 功能改进 2: Summary 图改为叠加显示

**改进背景**：
用户反馈希望在 Summary 图中看到所有 rod 的直方图叠加在同一张图上，而不是网格布局，以便更直观地比较不同 rod 的性能差异。

**实现方案**：
1. **修改 `plot_summary_grid()` 函数**：
   - 从网格布局（GridSpec）改为单图叠加
   - 使用 matplotlib 的 tab20 色图为每个 rod 分配颜色
   - 所有直方图曲线绘制在同一坐标系
   - 峰位置用 X 标记
   - 图例显示每个 rod 的 (iy, iz) 和峰数

2. **颜色映射**：
   - 使用 `plt.colormaps.get_cmap('tab20')` 替代已弃用的 `cm.get_cmap()`
   - 自动根据 rod 数量重采样颜色
   - 每个 rod 获得唯一颜色

3. **图例布局**：
   - 放置在图表右侧（`bbox_to_anchor=(1, 0.5)`）
   - 显示格式：`Rod(iy=X,iz=Y) P=峰数`
   - 半透明背景便于阅读

**视觉效果对比**：
- **旧版（网格）**：每个 rod 独立子图，难以比较
- **新版（叠加）**：所有 rod 在同一坐标系，峰位置对比一目了然

**代码示例**：
```python
# 为每个 rod 绘制叠加曲线
for idx, res in enumerate(results):
    color = colors(idx)
    label = f"Rod(iy={iy},iz={iz}) P={res['n_peaks']}"
    ax.plot(hist_x, hist_y_smooth, linewidth=1.5, alpha=0.7, 
           color=color, label=label)
    ax.scatter(hist_x[peaks], hist_y_smooth[peaks], s=30, 
              color=color, marker='x')

# 图例放在右侧
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
```

### Bug 修复 3: Matplotlib 弃用警告

**问题**：使用了已弃用的 `cm.get_cmap()` 函数
**修复**：改用 `plt.colormaps.get_cmap('tab20').resampled(n)`

### Bug 修复 4: 变量作用域错误

**问题**：在函数内部错误位置导入 plt，导致 UnboundLocalError
**修复**：移除函数内部的 import 语句，使用全局导入的 matplotlib.pyplot

## 实际运行测试结果

### 测试环境
- Python 3.14
- Windows 10
- 依赖库：pandas 3.0.0, numpy 2.4.1, matplotlib 3.10.8, scipy 1.17.0

### 测试数据
- 配置数量：3 个（001_Ny1, 002_Ny2, 003_Ny3）
- 晶体阵列：1×3, 2×3, 3×3
- 总 rod 数：18 个
- 每个配置 8 个线程文件

### 测试结果
```
[Processing] 001_Ny1
  [Merge] Created AnaEx01_nt_PhotonLRPerRod.csv (300 rows)
  3 rods analyzed: 10, 10, 10 peaks

[Processing] 002_Ny2  
  [Merge] Created AnaEx01_nt_PhotonLRPerRod.csv (600 rows)
  6 rods analyzed: 5-10 peaks

[Processing] 003_Ny3
  [Merge] Created AnaEx01_nt_PhotonLRPerRod.csv (900 rows)
  9 rods analyzed: 5-11 peaks

Total rods analyzed: 18
Execution time: ~20 seconds
```

### 生成文件
- 18 个单独 rod 直方图 PNG
- 3 个配置 summary 叠加图
- 1 个全局汇总 CSV
- 所有文件符合预期

## 后续优化方向

可能的改进点：
1. 添加多进程并行处理加速批量分析
2. 支持更多统计指标（Crosstalk、能量分辨率等）
3. 交互式 HTML 报告生成
4. 参数优化建议功能（基于 PVR 和 FWHM）
5. 支持自定义图表样式和颜色主题
6. 添加归一化选项（便于不同事件数的 rod 对比）

## 参考资料

- scipy.signal.find_peaks 文档
- scipy.optimize.curve_fit 文档
- matplotlib colormaps 文档
- pandas CSV 处理最佳实践
- Geant4 CSV 输出格式规范

---

**影响范围**: 新增数据分析工具，不影响现有 Geant4 仿真流程

**测试状态**: ✅ 实际数据测试通过，功能完整

**最终版本**: v1.1 (2026-01-31)

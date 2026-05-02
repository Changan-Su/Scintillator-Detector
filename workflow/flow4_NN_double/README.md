# flow4_NN_double

**双点光源**分支的批量仿真脚本。每个 config 包含 **两个源位置** A 和 B，
按 1:2 ~ 2:1 之间的随机比例瓜分总 `beamOn`（默认 35000）。
**单进程一次 `/run/beamOn`** 同时覆盖两点（C++ 端按 `fraction_a` 概率分流），
不需要事后合并。

## 与 flow2 / flow3_NN 的关系

| 阶段 | flow2 / flow3_NN（单点） | flow4_NN_double（双点） |
|------|--------------------------|--------------------------|
| 仿真 | `run_batch_sigma_position.py` 一个 config 一次 beamOn | `run_batch_double_point.py` 一个 config **同样**一次 beamOn，但 C++ 内部按概率从 A 或 B 发射 |
| 合并 | `split_events.py` 切大批量 | 不需要 |
| 重建 | `Histo10_Cubic.py` | 同左，直接对 `Results/<config>/` 跑即可（注意单点算法对双点意义有限） |

## 单位约定

**所有坐标和距离都是 mm**（Geant4 内部默认长度单位）。
`/source/fp_source` 和 `/source/fp_source_b` 接收的是裸 double，按 mm 解释。
脚本里的 `--random-min/--max`、`--min-separation`、`--max-separation`、`--point-a/-b`
也都按 mm 传入；`metadata.csv` 字段后缀就是 `_mm`。

## C++ 改动（已并入此分支）

`PrimaryGeneratorAction` 新增两个 messenger 命令：

| 命令 | 含义 | 默认 |
|------|------|------|
| `/source/fp_source_b X Y Z` | 第二个点位置（mm） | `(0, 0, 0)` |
| `/source/fraction_a F`      | 每个 primary 从 A 点发射的概率（B 点概率 = `1-F`） | `1.0` |

`fraction_a == 1.0` 时退化为单点行为，所有 flow1/flow2/flow3 旧 macro 完全不受影响。
`HistoManager` 也始终把 `source.fp_source_b_*_mm` 和 `source.fraction_a` 写进
`metadata.csv`，单点跑出来的 `fraction_a` 就是 1.0。

合并后的 metadata.csv 会包含所有理论真值：
```
source.fp_source_x_mm, source.fp_source_y_mm, source.fp_source_z_mm    ← A 点理论位置
source.fp_source_b_x_mm, source.fp_source_b_y_mm, source.fp_source_b_z_mm  ← B 点理论位置
source.fraction_a                                                       ← A 占比
```

## 用法

### 1. 检查计划（不真跑仿真）

```powershell
uv run python workflow/flow4_NN_double/run_batch_double_point.py `
    --random-count 5 --random-seed 42 --dry-run
```

### 2. 默认随机（10 个 config，±10 mm 立方体内）

```powershell
uv run python workflow/flow4_NN_double/run_batch_double_point.py `
    --random-count 10 --random-seed 42
```

### 3. 限制两点距离

```powershell
uv run python workflow/flow4_NN_double/run_batch_double_point.py `
    --random-count 100 `
    --min-separation 2 --max-separation 8 `
    --random-seed 1234
```

### 4. 显式两点（手动指定位置）

```powershell
uv run python workflow/flow4_NN_double/run_batch_double_point.py `
    --mode explicit --point-a 3 -2 5 --point-b -4 6 -1 `
    --fraction-a 0.4
```

不传 `--fraction-a` 时也会从 `[--ratio-min, --ratio-max]` 中随机一个。

## 关键参数

### 通用

| 参数 | 说明 | 默认 |
|------|------|------|
| `--mode` | `random` 或 `explicit` | `random` |
| `--beam-on-total N` | 单 config 总 beamOn | `35000` |
| `--ratio-min/--ratio-max` | A 占总 beamOn 的比例区间（1/3 ↔ 2/3 即 1:2 ~ 2:1） | `0.333 / 0.667` |
| `--sigma` | surfaceSigma 固定值 | `0.3` |
| `--results-dir` | Geant4 输出根 | `Results` |
| `--summary-csv` | 真值汇总 CSV 路径 | `<results-dir>/double_point_ground_truth.csv` |
| `--prefix-base` | config 名前缀 | `DP` |
| `--dry-run` | 只打印计划 | `False` |

### random 模式

| 参数 | 说明 | 默认 |
|------|------|------|
| `--random-count N` | 双点 config 总数 | `10` |
| `--random-min/--random-max` | xyz 立方体采样范围（mm） | `-10.0 / 10.0` |
| `--min-separation` | A、B 之间最小 3D 距离（mm），0=不限制 | `0` |
| `--max-separation` | A、B 之间最大 3D 距离（mm），0=不限制 | `0` |
| `--random-precision` | 坐标小数位 | `4` |
| `--random-seed` | RNG 种子 | `None` |

### explicit 模式

| 参数 | 说明 |
|------|------|
| `--point-a X Y Z` | A 点坐标（mm，必填） |
| `--point-b X Y Z` | B 点坐标（mm，必填） |
| `--fraction-a F`  | A 点比例（可选，省略则随机） |

## 产物

```
Results/
  DP_0001_S0p3_A2p79_m9p5_m4p5_Bm5p54_4p73_3p53_<ts>/
    AnaEx01_nt_PhotonFaceBlockEvent_t*.csv
    metadata.csv                   # 含 source.fp_source_b_* 和 source.fraction_a
  DP_0002_.../ ...
  double_point_ground_truth.csv    # 一行一个 config 的真值汇总
```

`double_point_ground_truth.csv` 列：

```
config_name, ax_mm, ay_mm, az_mm, bx_mm, by_mm, bz_mm,
fraction_a, expected_n_a, expected_n_b, separation_mm,
sigma, beam_on_total, results_subdir
```

注：`expected_n_a/b` 是按 fraction_a 算出来的期望值，**实际**每个事件是
C++ 端 `G4UniformRand() < fraction_a` 决定的，所以真实数量会有 √N 量级
的二项分布涨落（35000 × 0.5 ± 90 左右）。

## 可视化训练数据分布

跑完一个 flow4 扫描后，可以用 `visualize_points.py` 把所有 (A, B) 真值
位置画成 3D 立方体散点图（自动检测 `DP_*/metadata.csv`，
没有就回退读 `double_point_ground_truth.csv`）：

```powershell
uv run python workflow/flow4_NN_double/visualize_points.py Results_flow4_trainingdata --no-show --projections
```

产物在文件夹根下：`double_point_visualization.png`（3D 总览，A 红 B 蓝、
带 A→B 连线、25 mm 晶体线框），加 `--projections` 还会多出
`*_xy.png / *_xz.png / *_yz.png` 三张 2D 投影。常用开关：
`--alpha 0.4` 降透明度看密集分布、`--no-lines` 只看散点、
`--cube-size 20` 手动改立方体边长。

## 下游

```powershell
uv run python workflow/flow3_NN/Histo10_Cubic.py Results
```

经典单点重建算法（half_side / linear / mlp_*）默认会把双点拟合到加权重心
附近，并不是真正的双点定位。真正的双点重建需要后续训练专用 NN
（待加 `build_dataset_double.py` 等）。

## Qt 可视化（看双点光源效果）

项目根有 `vis_double_point.mac`（也在 `build/Release/` 下放了一份），
里面已经把双点参数和坐标标签都设好。

**1. 启动 Qt 交互模式** — 不带任何参数运行 exe：

```powershell
cd D:\Geant4\Projects\Scintillator-Detector-Continuous\build\Release
.\exampleB1.exe
```

启动后 Geant4 自动跑 `init_vis.mac`（`geometry.mac` → `/run/initialize` → `vis.mac`），
晶体几何出现在 Qt OpenGL 窗口里。

**2. 在 Qt 底部 Session 输入框里执行双点 macro**：

```
/control/execute vis_double_point.mac
```

里面的关键命令：

```
/source/mode optical
/source/distribution Point
/source/fp_source     3 -2  5     # A 点 (mm)
/source/fp_source_b  -4  6 -1     # B 点 (mm)
/source/fraction_a    0.4         # A 概率（B = 0.6）
/vis/scene/add/text  3 -2  5 mm 18 4 4 A
/vis/scene/add/text -4  6 -1 mm 18 4 4 B
/run/beamOn 200
```

会看到两簇黄色 optical photon（vis.mac 里 `opticalphoton yellow` 已设）从两个不同位置射出，
旁边带有 A/B 文字标记。

**3. 调参数 / 反复跑** — 直接在 Qt Session 里改单条命令再 `/run/beamOn`：

```
/source/fp_source     0  0  0
/source/fp_source_b  10  0  0
/source/fraction_a   0.5
/run/beamOn 200
```

清掉旧轨迹：

```
/vis/viewer/clearTransientStore
```

**事件数建议**：optical 模式下每 event = 1 photon，画 trajectory 比较吃内存。看相对位置和分布
200–1000 足够；数万就会卡。

> ⚠️ 必须用 **本分支重新编译过的** `exampleB1.exe`。旧 exe 不认识 `/source/fp_source_b`
> 和 `/source/fraction_a`，会报 "command not found"。

---

## 设计说明：为什么单进程？

`HistoManager::Book()` 在 `BeginOfRunAction` 打开 CSV，
`EndOfRunAction` 调用 `CloseFile()`。如果同一进程里连续两次
`/run/beamOn`，第二次会**覆盖**第一次的 CSV。

之前考虑过两次 subprocess（A、B 各跑一次再 Python 合并），但用户要求
"必须在同一个模拟里面跑"。所以改在 C++ 里加 `fp_source_b` + `fraction_a`，
让单次 beamOn 内部就把两点同时跑掉，CSV 输出格式和单点完全一致。

---

## NN 双点重建（训练 + 推理）

5 个脚本组成一条 build → train → eval → predict 的管线，全部使用 mm。
loss 对 A↔B 排列不变（`min(identity, swap)`），训练自动学到的 (A, B) 顺序
任意；评估和推理时按"最优匹配"对齐到真值。

```powershell
# 1. 构建数据集（可同时合并多个 separation regime）
uv run python workflow/flow4_NN_double/build_dataset_double.py `
    --input-roots Results_flow4_trainingdata `
                  Results_flow4_trainingdata_sep5mm `
                  Results_flow4_trainingdata_sep3mm `
    --photons-per-sample 2000 --normalize-counts

# 2. 训练（GPU 大 batch、AdamW、ReduceLROnPlateau、early stop）
uv run python workflow/flow4_NN_double/train_double.py --epochs 300

# 3. 测试集评估（permutation-aware MAE / RMSE / separation 误差）
uv run python workflow/flow4_NN_double/evaluate_double.py

# 4. 单 config 推理（输出 nn_double_predicted.csv 到该 config 目录）
uv run python workflow/flow4_NN_double/predict_double.py `
    --config-dir Results_flow4_trainingdata/DP_0001_... --normalize-counts
```

产物在 `workflow/flow4_NN_double/artifacts/`：
`dataset.npz` / `best.pt` / `norm.npz` / `split_idx.npz` /
`loss_curve.png` / `test_metrics.csv` / `test_scatter.png`。

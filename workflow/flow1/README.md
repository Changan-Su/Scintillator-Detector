# flow1

`flow1` 是当前主分析流程的一份独立副本，针对六面 SiPM 的这条链路：

1. 批量跑仿真到 `Results/`
2. 从 `Results/` 生成 6 面光子计数 heatmap 和重建点位
3. 对比模拟真值点位和重建点位

本目录内包含的脚本：

- `run_batch_sigma_position.py`
- `Histo10_Cubic.py`
- `analyze_position_accuracy.py`

## 脚本使用顺序

按下面顺序执行：

1. `python workflow\flow1\run_batch_sigma_position.py`
2. `python workflow\flow1\Histo10_Cubic.py Results`
3. `python workflow\flow1\analyze_position_accuracy.py Output\SiPM6_Output_<timestamp>`

## 每个脚本的用途

### `run_batch_sigma_position.py`

用途：

- 批量调用 `build\Release\exampleB1.exe`
- 扫描 `/detector/surfaceSigma`
- 扫描 `/source/fp_source x y z`
- 为每次运行自动设置 `/results/prefix`
- 生成 `Results/<config>/...` 原始仿真输出

典型输入：

- `build\Release\exampleB1.exe`
- `geometry.mac`

典型输出：

- `Results/<config>/metadata.csv`
- `Results/<config>/AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`

### `Histo10_Cubic.py`

用途：

- 读取每个 `Results/<config>/` 下的 `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`
- 自动合并线程 CSV
- 聚合为六面 4x4 SiPM 计数矩阵
- 输出热图和多种重建算法结果

输入：

- `Results/<config>/AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`
- 可选 `Results/<config>/metadata.csv`

输出目录：

- `Output/SiPM6_Output_<timestamp>/<config>/`

每个配置目录下会生成：

- `merged_event.csv`
- `merged_face_jk.csv`
- `reconstructed_position.csv`
- `metadata.csv`
- `sipm_6faces_heatmap.png`

### `analyze_position_accuracy.py`

用途：

- 读取 `metadata.csv` 中的模拟真值点位
- 读取 `reconstructed_position.csv` 中的重建点位
- 按算法统计误差
- 生成对比散点图和精度汇总表

输入：

- `Output/.../<config>/metadata.csv`
- `Output/.../<config>/reconstructed_position.csv`

输出：

- `accuracy_summary.csv`
- `scatter_<algorithm>.png`

## 最短执行示例

```powershell
cmake --build build --config Release
python workflow\flow1\run_batch_sigma_position.py --beam-on 5000
python workflow\flow1\Histo10_Cubic.py Results
python workflow\flow1\analyze_position_accuracy.py Output\SiPM6_Output_<timestamp>
```

把 `<timestamp>` 替换成 `Histo10_Cubic.py` 实际创建出的输出目录名。

## 数据流

```text
workflow\flow1\run_batch_sigma_position.py
  -> Results/<config>/
     -> metadata.csv
     -> AnaEx01_nt_PhotonFaceBlockEvent_t*.csv

workflow\flow1\Histo10_Cubic.py Results
  -> Output/SiPM6_Output_<timestamp>/<config>/
     -> merged_event.csv
     -> merged_face_jk.csv
     -> reconstructed_position.csv
     -> metadata.csv
     -> sipm_6faces_heatmap.png

workflow\flow1\analyze_position_accuracy.py Output/SiPM6_Output_<timestamp>
  -> Output/SiPM6_Output_<timestamp>/accuracy_summary.csv
  -> Output/SiPM6_Output_<timestamp>/scatter_<algorithm>.png
```

## 说明

- 这套流程使用的是 `PhotonFaceBlockEvent` 六面数据。
- `merge_photon_lr_perrod.py` 和 `Python_Scripts/Histo9.py` 不属于这条主流程，它们是旧的 `PhotonLRPerRod` 分析链路。
- 这里复制的是当前仓库里的工作副本；如果根目录脚本后续更新，`flow1` 里的副本不会自动同步。

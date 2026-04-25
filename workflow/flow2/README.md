# flow2

`flow2` 在 flow1 的基础上增加了**事件拆分**步骤：
对一次大批量仿真（如 10 万个事件）的原始数据，
按可配置的事件数切分为若干独立批次（如 10 组 × 1 万），
然后对每个批次分别重建，用于观察统计误差分布。

本目录内包含的脚本：

- `run_batch_sigma_position.py`  （与 flow1 相同）
- `split_events.py`              （新增：事件拆分）
- `Histo10_Cubic.py`             （与 flow1 相同）
- `analyze_position_accuracy.py` （与 flow1 相同）

## 脚本使用顺序

```powershell
python workflow\flow2\run_batch_sigma_position.py --beam-on 100000
python workflow\flow2\split_events.py Results --batch-size 10000
python workflow\flow2\Histo10_Cubic.py Results_split
python workflow\flow2\analyze_position_accuracy.py Output\SiPM6_Output_<timestamp>
```

把 `<timestamp>` 替换成 `Histo10_Cubic.py` 实际创建出的输出目录名。

## 每个脚本的用途

### `run_batch_sigma_position.py`

与 flow1 完全相同。用于批量调用仿真并将原始数据写入 `Results/<config>/`。

典型输入：
- `build\Release\exampleB1.exe`
- `geometry.mac`

典型输出：
- `Results/<config>/metadata.csv`
- `Results/<config>/AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`

---

### `split_events.py`（新增）

用途：

- 读取 `Results/<config>/` 下的所有线程 CSV
- 合并后按 EventID 顺序切分为若干等大批次
- 每批次独立写到 `Results_split/<config>_batch_<N>/`，并复制 `metadata.csv`

参数：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `results` | Results 根目录 | （必填） |
| `--batch-size N` | 每批次的事件数 | `10000` |
| `--output PATH` | 拆分输出根目录 | `Results_split` |
| `--dry-run` | 只打印计划，不写文件 | `False` |

示例：

```powershell
# 查看拆分计划（不写文件）
python workflow\flow2\split_events.py Results --batch-size 10000 --dry-run

# 实际拆分
python workflow\flow2\split_events.py Results --batch-size 10000

# 自定义批次大小和输出目录
python workflow\flow2\split_events.py Results --batch-size 5000 --output Results_5k
```

输入目录结构：

```
Results/
  S0p3_X0_Y0_Z0/
    AnaEx01_nt_PhotonFaceBlockEvent_t0.csv
    AnaEx01_nt_PhotonFaceBlockEvent_t1.csv
    metadata.csv
```

输出目录结构（`--batch-size 10000`，共 10 万事件）：

```
Results_split/
  S0p3_X0_Y0_Z0_batch_0001/
    AnaEx01_nt_PhotonFaceBlockEvent.csv   (事件 1..10000)
    metadata.csv
  S0p3_X0_Y0_Z0_batch_0002/
    AnaEx01_nt_PhotonFaceBlockEvent.csv   (事件 10001..20000)
    metadata.csv
  ...
  S0p3_X0_Y0_Z0_batch_0010/
    AnaEx01_nt_PhotonFaceBlockEvent.csv   (事件 90001..100000)
    metadata.csv
```

---

### `Histo10_Cubic.py`

与 flow1 完全相同。输入改为 `Results_split` 而非 `Results`，
脚本会遍历所有批次子目录并分别重建。

```powershell
python workflow\flow2\Histo10_Cubic.py Results_split
```

每个批次子目录对应 `Output/SiPM6_Output_<timestamp>/` 下的同名子目录，
各自生成 `reconstructed_position.csv`。

---

### `analyze_position_accuracy.py`

与 flow1 完全相同。汇总所有批次的重建误差，
额外体现批次间的统计分散（std 可反映不同事件数下的统计涨落）。

```powershell
python workflow\flow2\analyze_position_accuracy.py Output\SiPM6_Output_<timestamp>
```

## 数据流

```text
workflow\flow2\run_batch_sigma_position.py --beam-on 100000
  -> Results/<config>/
       metadata.csv
       AnaEx01_nt_PhotonFaceBlockEvent_t*.csv   (100000 events)

workflow\flow2\split_events.py Results --batch-size 10000
  -> Results_split/<config>_batch_0001/
       AnaEx01_nt_PhotonFaceBlockEvent.csv      (10000 events)
       metadata.csv
  -> Results_split/<config>_batch_0002/ ...
  -> ...（共 10 批）

workflow\flow2\Histo10_Cubic.py Results_split
  -> Output/SiPM6_Output_<timestamp>/<config>_batch_0001/
       merged_event.csv
       merged_face_jk.csv
       reconstructed_position.csv
       metadata.csv
       sipm_6faces_heatmap.png
  -> Output/SiPM6_Output_<timestamp>/<config>_batch_0002/ ...

workflow\flow2\analyze_position_accuracy.py Output/SiPM6_Output_<timestamp>
  -> Output/SiPM6_Output_<timestamp>/accuracy_summary.csv
  -> Output/SiPM6_Output_<timestamp>/scatter_<algorithm>.png
```

## 与 flow1 的区别

| | flow1 | flow2 |
|-|-------|-------|
| 仿真次数 | 每个参数组合独立跑 | 每个参数组合跑一次（大批量） |
| 重建次数 | 每次仿真 1 次重建 | 每次仿真拆为 N 批，各重建 1 次 |
| 新增脚本 | — | `split_events.py` |
| `Histo10_Cubic.py` 输入 | `Results/` | `Results_split/` |
| 典型应用 | 扫描参数空间 | 评估统计误差 |

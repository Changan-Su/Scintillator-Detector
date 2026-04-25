# Scintillator-Detector-Continuous

Geant4-based scintillator detector simulation with a Python post-processing pipeline for:

- sweeping source position and optical surface sigma
- collecting per-run CSV output under `Results/`
- generating 6-face SiPM photon-count heatmaps
- reconstructing source position from the 6-face response
- comparing reconstructed positions against the true simulated positions

The main executable target is `exampleB1`.

## Build

### Windows

```powershell
cmake -S . -B build -DGeant4_DIR=D:/Geant4/geant4-install/lib/cmake/Geant4 -DCMAKE_PREFIX_PATH=D:/Qt2/5.15.2/msvc2019_64/lib/cmake
cmake --build build --config Release
```

### Run one macro

```powershell
build\Release\exampleB1.exe run4.mac
```

Useful alternatives:

- `run3.mac`: faster smoke test
- `run4.mac`: typical batch-style run
- `run_vis.bat`: visualization helper

## Current Python Workflow

This is the current six-face SiPM workflow. If your goal is:

- start from many `Results/<config>/` folders
- make per-config heatmaps
- reconstruct positions
- compare true source positions with reconstructed positions

then the script order is:

1. `run_batch_sigma_position.py`
2. `Histo10_Cubic.py`
3. `analyze_position_accuracy.py`

### Step 1: batch-run simulations into `Results/`

```powershell
python run_batch_sigma_position.py
```

What it does:

- runs `build\Release\exampleB1.exe` repeatedly
- sweeps `/detector/surfaceSigma`
- sweeps `/source/fp_source x y z`
- sets `/results/prefix` automatically for each run
- creates one subdirectory per run under `Results/`

Typical per-run output under `Results/<config>/`:

- `metadata.csv`
- `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`
- other Geant4 CSV outputs

Useful options:

```powershell
python run_batch_sigma_position.py --dry-run
python run_batch_sigma_position.py --beam-on 50000
python run_batch_sigma_position.py --sigma-start 0.3 --sigma-end 0.7 --sigma-step 0.1
python run_batch_sigma_position.py --x-start -12.5 --x-end 12.5 --x-step 6.25
```

### Step 2: merge `Results/` and generate heatmaps + reconstructed positions

```powershell
python Histo10_Cubic.py Results
```

What it does for each `Results/<config>/` folder:

- merges `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`
- writes `merged_event.csv`
- aggregates photon counts by `(Face, j, k)` into `merged_face_jk.csv`
- generates `sipm_6faces_heatmap.png`
- computes several reconstruction algorithms and writes `reconstructed_position.csv`
- copies `metadata.csv` into the output folder

Default output location:

```text
Output/SiPM6_Output_<timestamp>/<config>/
```

Files produced per config:

- `merged_event.csv`
- `merged_face_jk.csv`
- `reconstructed_position.csv`
- `metadata.csv`
- `sipm_6faces_heatmap.png`

Useful options:

```powershell
python Histo10_Cubic.py Results --output Output\MyBatch
python Histo10_Cubic.py Results --cmap plasma
python Histo10_Cubic.py Results --vmin 0 --vmax 2000
```

Important note:

- `Histo10_Cubic.py` already merges the thread CSVs by itself.
- You do not need to run `merge_photon_lr_perrod.py` before `Histo10_Cubic.py`.

### Step 3: compare true source position vs reconstructed position

```powershell
python analyze_position_accuracy.py Output\SiPM6_Output_<timestamp>
```

What it reads from each output subfolder:

- `metadata.csv`
- `reconstructed_position.csv`

What it produces:

- `accuracy_summary.csv`
- `scatter_<algorithm>.png`

What it computes:

- true position from `metadata.csv` using `source.fp_source_x_mm`, `source.fp_source_y_mm`, `source.fp_source_z_mm`
- reconstructed position from each row in `reconstructed_position.csv`
- per-algorithm `dx`, `dy`, `dz`, and 3D distance error `d`
- per-algorithm bias and standard deviation summary

Useful options:

```powershell
python analyze_position_accuracy.py Output\SiPM6_Output_20260314_220905
python analyze_position_accuracy.py Output\SiPM6_Output_20260314_220905 --no-plot
python analyze_position_accuracy.py Output\SiPM6_Output_20260314_220905 --output Output\SiPM6_Output_20260314_220905\accuracy_summary_custom.csv
```

## What Each Current Script Is For

### `run_batch_sigma_position.py`

Purpose:

- batch driver for simulation
- parameter sweep over `surfaceSigma` and source point `(x, y, z)`

Input:

- built executable `build\Release\exampleB1.exe`
- `geometry.mac`

Output:

- many run folders under `Results/`

Use it when:

- you need fresh simulation data

### `Histo10_Cubic.py`

Purpose:

- main six-face analysis script
- converts raw `PhotonFaceBlockEvent` CSVs into heatmaps and reconstructed positions

Input:

- `Results/<config>/AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`
- optional `Results/<config>/metadata.csv`

Output:

- `Output/SiPM6_Output_<timestamp>/<config>/...`

Use it when:

- you want the per-config 6-face photon heatmap
- you want reconstructed positions for each config

### `analyze_position_accuracy.py`

Purpose:

- batch accuracy evaluation across all processed configs

Input:

- `Output/.../<config>/metadata.csv`
- `Output/.../<config>/reconstructed_position.csv`

Output:

- `accuracy_summary.csv`
- `scatter_<algorithm>.png`

Use it when:

- you want to compare simulated source positions with reconstructed positions

### `merge_photon_lr_perrod.py`

Purpose:

- merges old `AnaEx01_nt_PhotonLRPerRod_t*.csv` files into one CSV

Input:

- `AnaEx01_nt_PhotonLRPerRod_t*.csv`

Output:

- merged `AnaEx01_nt_PhotonLRPerRod_merged.csv`

Use it when:

- you are still using the older left/right per-rod analysis pipeline
- you are not using the current six-face `PhotonFaceBlockEvent` workflow

### `Python_Scripts/Histo9.py`

Purpose:

- legacy DOI analysis for `PhotonLRPerRod`
- analyzes per-rod left/right photon counts and detects DOI peaks

Use it when:

- you want DOI-style rod-by-rod analysis from the older `PhotonLRPerRod` data

Not part of the current six-face heatmap + reconstruction workflow.

## Current Data Flow

```text
run_batch_sigma_position.py
  -> Results/<config>/
     -> metadata.csv
     -> AnaEx01_nt_PhotonFaceBlockEvent_t*.csv

Histo10_Cubic.py Results
  -> Output/SiPM6_Output_<timestamp>/<config>/
     -> merged_event.csv
     -> merged_face_jk.csv
     -> reconstructed_position.csv
     -> metadata.csv
     -> sipm_6faces_heatmap.png

analyze_position_accuracy.py Output/SiPM6_Output_<timestamp>
  -> Output/SiPM6_Output_<timestamp>/accuracy_summary.csv
  -> Output/SiPM6_Output_<timestamp>/scatter_<algorithm>.png
```

## Shortest End-to-End Example

```powershell
cmake --build build --config Release
python run_batch_sigma_position.py --beam-on 5000
python Histo10_Cubic.py Results
python analyze_position_accuracy.py Output\SiPM6_Output_<timestamp>
```

Replace `<timestamp>` with the actual output directory created by `Histo10_Cubic.py`.

## Notes

- The current six-face workflow is based on `PhotonFaceBlockEvent`, not `PhotonLRPerRod`.
- `metadata.csv` is the source of truth for the simulated source position and detector settings for each run.
- `reconstructed_position.csv` stores multiple algorithms. `analyze_position_accuracy.py` evaluates all of them automatically.
- `sipm_6faces_heatmap.png` is the per-config photon-count heatmap. The true-vs-reconstructed comparison is produced separately by `scatter_<algorithm>.png`.

## Repository Layout

```text
src/                    C++ source files
include/                C++ headers
build/                  build output
Results/                raw Geant4 run outputs
Output/                 processed Python analysis outputs
run_batch_sigma_position.py
Histo10_Cubic.py
analyze_position_accuracy.py
merge_photon_lr_perrod.py
Python_Scripts/         older analysis scripts
```

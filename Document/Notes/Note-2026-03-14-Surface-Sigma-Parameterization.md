# Note 2026-03-14：Surface Sigma 参数化与批量扫参

本文整理「将晶体表面 Sigma 从硬编码改为宏可调、批处理可循环」的实现方案与知识点，便于日后复现或扩展类似参数。

---

## 1. 需求背景

- **问题**：`DetectorConstruction.cc` 中 `Surface_Sigma` 原本是局部硬编码（如 `0.7`），每次改值都要改源码、重新编译，无法批量扫参。
- **目标**：通过宏命令 `/detector/surfaceSigma` 在运行时设定，配合 `geometry.mac` 和 `run_batch.bat` 实现全自动按设定数值循环跑结果。

---

## 2. 方案架构

```
geometry.mac (/detector/surfaceSigma 0.7)
        ↓
DetectorMessenger::SetNewValue()  →  fDetector->SetSurfaceSigma(v)
        ↓
DetectorConstruction::fSurfaceSigma
        ↓
Construct() 中 crystalsurface->SetSigmaAlpha(fSurfaceSigma)
```

- **数据流**：宏命令 → Messenger 解析 → DetectorConstruction 成员变量 → 几何构建时使用。
- **时机**：`geometry.mac` 在 `run4.mac` 中于 `/run/initialize` 之前执行，因此 `Construct()` 被调用时已拿到正确值。

---

## 3. 实现步骤拆解

### 3.1 `include/DetectorConstruction.hh`

| 步骤 | 内容 | 说明 |
|------|------|------|
| 1 | 增加 `G4double GetSurfaceSigma() const;` | 只读访问，供调试或后续分析 |
| 2 | 增加 `void SetSurfaceSigma(G4double v);` | 参数校验：`v ∈ [0,1]` 否则用默认值 |
| 3 | 增加 `G4double fSurfaceSigma = 0.5;` | 私有成员，默认值需与 `SetSurfaceSigma` 的 fallback 一致 |

**注意**：Geant4 标量类型为 `G4double`（小写 d），不是 `G4Double`。写错会导致编译失败。

### 3.2 `include/DetectorMessenger.hh`

| 步骤 | 内容 | 说明 |
|------|------|------|
| 1 | 增加 `G4UIcmdWithADouble* fSurfaceSigmaCmd = nullptr;` | 与其他 `*Cmd` 成员声明同列 |

### 3.3 `src/DetectorMessenger.cc`

| 步骤 | 内容 | 说明 |
|------|------|------|
| 1 | 构造函数中创建命令 | `fSurfaceSigmaCmd = new G4UIcmdWithADouble("/detector/surfaceSigma", this);` |
| 2 | 设置参数名、范围、默认值 | `SetParameterName`, `SetRange("Sigma >= 0. && Sigma <= 1.")`, `SetDefaultValue` |
| 3 | 析构函数中 `delete fSurfaceSigmaCmd;` | 防止内存泄漏 |
| 4 | `SetNewValue()` 中增加分支 | `else if (command == fSurfaceSigmaCmd) fDetector->SetSurfaceSigma(...)` |

**顺序**：`SetNewValue` 中 `fSurfaceSigmaCmd` 的 `else if` 应放在 `fUpdateCmd` 之前，否则 `fUpdateCmd` 会先匹配到。

### 3.4 `src/DetectorConstruction.cc`

| 步骤 | 内容 | 说明 |
|------|------|------|
| 1 | 将 `G4double Surface_Sigma = 0.7;` 改为 `G4double Surface_Sigma = fSurfaceSigma;` | 或直接 `crystalsurface->SetSigmaAlpha(fSurfaceSigma)` |

---

## 4. 物理与 Geant4 知识点

### 4.1 SetSigmaAlpha 含义

- **用途**：`G4OpticalSurface::SetSigmaAlpha(sigma)` 定义光学表面的**粗糙度**（sigma alpha）。
- **取值**：通常 0～1，0 表示完全光滑，1 表示完全粗糙。
- **影响**：`ground` 表面下，sigma 控制反射/透射的角分布，影响光收集效率。

### 4.2 光学表面类型

- `SetType(dielectric_dielectric)`：介质-介质界面。
- `SetModel(unified)`：统一光学模型。
- `SetFinish(ground)`：磨砂表面，需配合 `SetSigmaAlpha` 使用。

### 4.3 UI 命令流程

1. 宏文件执行 `/detector/surfaceSigma 0.6`
2. Geant4 解析命令，找到 `DetectorMessenger` 注册的 `fSurfaceSigmaCmd`
3. 调用 `SetNewValue(command, "0.6")`
4. `GetNewDoubleValue("0.6")` 得到 `0.6`，传给 `SetSurfaceSigma`

---

## 5. 宏使用方式

在 `geometry.mac`（或批处理生成的 geometry 片段）中增加：

```
/detector/surfaceSigma 0.7
```

必须放在 `/run/initialize` 之前（`run4.mac` 已通过 `execute geometry.mac` 保证顺序）。

---

## 6. 与 run_batch 的集成（已实现）

`run_batch.bat` 与 `run_batch.sh` 已支持 Sigma 循环，配置项如下：

| 配置项 | 说明 | 示例 |
|--------|------|------|
| `LOOP_SURFACE_SIGMA` | 是否循环 Sigma | `true` / `false` |
| `SIGMA_START` / `SIGMA_END` / `SIGMA_STEP` | 扫描范围与步长 | `0.5` / `0.7` / `0.1` |
| `DEFAULT_SURFACE_SIGMA` | 不循环时使用的默认值 | `0.5` |
| `NAME_INCLUDE_SURFACE_SIGMA` | 是否在结果目录名中包含 Sigma | `true` / `false` |

- `generate_geometry_mac` 已增加第 11 个参数，输出 `/detector/surfaceSigma <值>`
- 文件夹命名：`_Sigma0p7` 表示 sigma=0.7（小数点用 `p` 替代，避免路径问题）

---

## 7. 常见坑：G4Double vs G4double

- **错误**：`void SetSurfaceSigma(G4Double v)`（大写 D）
- **原因**：Geant4 的 `globals.hh` 中只定义 `G4double`（小写），`G4Double` 未定义。
- **现象**：MSVC 报 `error C2061: 语法错误: 标识符"G4Double"`，以及后续 `v` 未声明。
- **正确**：`void SetSurfaceSigma(G4double v)`

---

## 8. 校验清单

- [ ] 编译通过：`cmake --build build --config Release --target exampleB1`
- [ ] 单次运行：`geometry.mac` 中写 `/detector/surfaceSigma 0.3`，运行 `run4.mac` 无报错
- [ ] 结果可区分：改不同 sigma 值，输出结果应有差异（如光子计数分布）

---

## 9. 相关文件

| 文件 | 作用 |
|------|------|
| `include/DetectorConstruction.hh` | 成员、setter/getter 声明 |
| `include/DetectorMessenger.hh` | 命令指针声明 |
| `src/DetectorConstruction.cc` | 使用 `fSurfaceSigma` 构建光学表面 |
| `src/DetectorMessenger.cc` | 命令创建、析构、`SetNewValue` 分发 |
| `geometry.mac` | 运行时传入 sigma 值 |
| `run_batch.bat` / `run_batch.sh` | 已支持 Sigma 循环 |

---

## 10. 常见误区与排查：为什么改了 Sigma 结果还一样

### 现象

- 在 `run5.mac` 中写了不同的 `/detector/surfaceSigma`，但输出热力图几乎不变。

### 根因

- `surfaceSigma` 在当前工程中用于几何构建阶段（`Construct()` 中 `SetSigmaAlpha`）。
- 若命令写在 `/run/initialize` 之后，几何已经构建完成，参数变化不会自动作用到当前几何。
- 仅改参数但不重建几何（未执行 `/detector/update`）时，本次 `beamOn` 仍使用旧几何。

### 正确用法

**方案 A（推荐）**：在初始化前设置 Sigma

```
/control/execute geometry.mac
/detector/surfaceSigma 0.1
/run/initialize
/run/beamOn 100000
```

**方案 B**：初始化后改值时，显式重建几何

```
/run/initialize
/detector/surfaceSigma 0.1
/detector/update
/run/beamOn 100000
```

### 快速检查清单

- [ ] `/detector/surfaceSigma` 是否在 `/run/initialize` 之前？
- [ ] 若在之后，是否执行了 `/detector/update`？
- [ ] `geometry.mac` 是否被当前 run 宏正确执行（`/control/execute geometry.mac`）？

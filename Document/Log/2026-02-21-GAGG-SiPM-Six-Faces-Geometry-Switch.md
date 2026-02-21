# Log: GAGG 唯一闪烁体、SiPM 六面 4×4、几何参数开关与可视化修复

**日期**: 2026-02-21  
**项目**: Scintillator-Detector-Continuous  

---

## 开发内容概览

1. **闪烁体仅用 GAGG**：移除 YSO，晶体与 Fillter 均使用 GAGG 材料及 GAGG-HL 光学/闪烁属性。  
2. **晶体阵列参数来源开关**：新增 bool 开关，可选“geometry.mac/UI 参数”或“手动参数”，并配套 `/detector/manual*` 等命令。  
3. **SiPM 六面 4×4 排布**：每个晶体 6 个面各贴 4×4 SiPM（像素 6 mm，间隙 0.2 mm），共 6×16=96 个 SiPM/晶体。  
4. **Qt 可视化六面显示修复**：通过唯一物理体积名与旋转矩阵常驻，使 6 个面的 SiPM 在 Qt/OpenGL 中均正确显示。  
5. **编译与编码**：关键注释改为英文、显式类型与括号匹配修正，消除 MSVC 代码页 936 下的解析错误。

---

## 实现方式

### 1) GAGG 唯一闪烁体

- **文件**: `src/DetectorConstruction.cc`
- 保留 GAGG 材料定义（Gd3Al2Ga3O12，密度 6.6 g/cm³）及 GAGG-HL 材料属性表（折射率 1.91，光产额 54 000 ph/MeV，衰减 150 ns，发射峰约 530 nm）。
- 删除原 YSO 用 MPT 块及 `YSO->SetMaterialPropertiesTable(...)`。
- `logicCrystal`、`logicFillter` 的 `G4LogicalVolume` 构造均改为使用 `GAGG` 材料。
- 晶体与 Fillter 的闪烁与光学行为统一由 GAGG 的 `gagg_mt` 决定。

### 2) 晶体阵列参数来源开关（geometry.mac vs 手动）

- **文件**: `include/DetectorConstruction.hh`
  - 新增 `fUseGeometryMac`（默认 true）及一套手动参数：`fManual_nx/ny/nz`、`fManual_crystal_gap`、`fManual_crystal_l/ly`、`fManual_fillter_ratio_*`、`fManual_fillter_pos_*`。
  - 新增对应 setter：`SetUseGeometryMac`、`SetManualArrayNx/Ny/Nz`、`SetManualCrystalGap/Size/SizeY`、`SetManualFillterRatio*`、`SetManualFillterPosRatio*`。
- **文件**: `include/DetectorMessenger.hh`、`src/DetectorMessenger.cc`
  - 新增命令：`/detector/useGeometryMac`（bool）、`/detector/manualArrayNx/Ny/Nz`、`/detector/manualCrystalGap`、`/detector/manualCrystalSize/SizeY`、`/detector/manualFillterRatioY/Z`、`/detector/manualFillterPosRatioY/Z`。
- **文件**: `src/DetectorConstruction.cc`
  - 在晶体阵列构建前，根据 `fUseGeometryMac` 选择使用“geometry.mac 参数”（`fPar_*`、`fCrystal_gap`、`fcrystal_l/ly`、Fillter 比例）或“手动参数”（`fManual_*`），再计算 `Crystal_nx/ny/nz`、`Crystal_gap`、`crystal_l/ly` 及 Fillter 比例，后续几何与 SiPM 放置均基于此组变量。

### 3) SiPM 六面 4×4 排布

- **文件**: `src/DetectorConstruction.cc`
- **单像素与间隙**：`sipm_l = 6 mm`，`SiPm_gap = 0.2 mm`，每面 4×4 像素，中心距 `pitch = 6.2 mm`（对应 25 mm 晶体单边 4 像素 + 5 间隙）。
- **几何**：SiPM 为薄盒，厚度 `sipm_t = 0.6 mm` 沿法向贴晶体表面，6×6 mm 在面内；`G4Box` 半边长顺序为 (厚度/2, 边长/2, 边长/2)，放置时通过旋转使厚度方向对准各面法向。
- **6 个面**：face 0~5 对应 +X、-X、+Y、-Y、+Z、-Z；每个面中心在晶体外表面外偏 `sipm_t/2`，面内用 `(j-1.5)*pitch`、`(k-1.5)*pitch` 做 4×4 网格偏移。
- **旋转**：+X 不旋转（nullptr）；-X 绕 Y 转 180°；±Y 绕 Z 转 ±90°；±Z 绕 Y 转 ±90°。旋转矩阵在循环外创建 5 个（`rotSiPM_mX/pY/mY/pZ/mZ`）并复用，**不在循环内 delete**，避免几何树中旋转失效导致仅一面显示。
- **Copy 号与命名**：`copyNoSiPM = crystalId*6*16 + face*16 + j*4 + k`；物理体积名 `SiPM_%d`（copyNoSiPM）保证每个 SiPM 唯一名，便于 Qt 场景树与可视化正确显示 6 面。

### 4) 编译与编码修复

- **文件**: `src/DetectorConstruction.cc`
- 将 GAGG、晶体阵列、Fillter、SiPM 及三重循环与 6 面 SiPM 段落中的**中文注释改为英文**，避免在 MSVC 代码页 936 下多字节字符被误解析为 `}` 等，导致“未声明的标识符”“非法 else”“缺少类型说明符”等连锁错误。
- `logicCrystal`、`solidCrystal` 改为显式类型 `G4LogicalVolume*`、`G4Box*`，避免在 YSO 移除后类型推断异常。
- 保持 `Construct()` 内大括号匹配正确，`fScoringVolume`、`fCrystal_*`、`flogicSiPM` 赋值与 `return physWorld` 均在函数内。

---

## 验证与注意事项

- **构建**：`cmake --build build --config Release` 通过，生成 `exampleB1.exe`；若出现 C4819，可将源文件存为“带 BOM 的 UTF-8”或保持当前英文注释。
- **运行**：geometry.mac 或 UI 设置晶体参数后，`/detector/update` 重建几何；手动模式时先 `/detector/useGeometryMac false`，再设 `/detector/manual*` 后 update。
- **可视化**：Qt 中 6 个面的 SiPM 应均可见；若曾出现“仅一面 16 个”，系旋转矩阵在循环内被 delete 导致，已通过“循环外创建、复用、不 delete”修复。

---

## 对既有开发内容的修正说明

- 此前 SiPM 仅在 ±X 两列放置且依赖已注释的 `sipm_l_ratio`，现改为每晶体 6 面 × 4×4，并统一 6 mm / 0.2 mm 规格。
- 晶体阵列参数曾仅来自 geometry.mac；现可通过 `fUseGeometryMac` 与 `/detector/manual*` 使用一套独立手动参数，便于脚本或固定配置测试。

---

## 涉及文件一览

| 文件 | 修改要点 |
|------|----------|
| `src/DetectorConstruction.cc` | GAGG 唯一、几何分支、SiPM 6 面 4×4、旋转常驻、注释/类型/括号 |
| `include/DetectorConstruction.hh` | `fUseGeometryMac`、`fManual_*` 及 setter |
| `src/DetectorMessenger.cc` | `/detector/useGeometryMac`、`/detector/manual*` 命令 |
| `include/DetectorMessenger.hh` | 上述命令成员指针 |

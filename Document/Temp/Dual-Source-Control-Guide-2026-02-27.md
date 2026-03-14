# Geant4 双发射源控制完整方案（基于当前项目）

## 1. 目标

在当前项目中实现并控制两把“枪”：

- 枪A：`gamma`（你现在已有）
- 枪B：`opticalphoton`（新增）

并通过宏命令控制模式：

- `gamma`：只发 gamma
- `optical`：只发 opticalphoton
- `both`：两者都发

> 关键概念：  
> 你不是要新增一个 Geant4 类，而是在同一个 `PrimaryGeneratorAction` 里新增一个 `G4ParticleGun*` 成员。

---

## 2. 当前项目现状（已对齐你的代码）

你当前代码结构是：

- `ActionInitialization` 只注册一个 `PrimaryGeneratorAction`
- `PrimaryGeneratorAction` 里目前只有 `fParticleGun`
- `run4mini.mac` 用 `/run/beamOn X`

这完全没问题。  
`/run/beamOn X` 只控制事件数，不控制“发哪把枪”。  
发哪把枪由 `PrimaryGeneratorAction::GeneratePrimaries()` 决定。

---

## 3. 实现方案总览

改两个文件：

1. `include/PrimaryGeneratorAction.hh`
2. `src/PrimaryGeneratorAction.cc`

核心动作：

1) 新增第二把枪 `fOpticalGun`  
2) 新增模式变量 `fSourceMode`  
3) 新增 messenger 命令 `/source/mode`  
4) 在 `GeneratePrimaries()` 里按模式调用对应 gun

---

## 4. 逐步修改（可直接照抄）

## Step A. 修改 `include/PrimaryGeneratorAction.hh`

### A1) 在 include 区新增

```cpp
#include "G4GenericMessenger.hh"
```

### A2) 在 `private:` 成员里新增/调整为

```cpp
G4ParticleGun* fParticleGun = nullptr;     // 原 gamma 枪
G4ParticleGun* fOpticalGun = nullptr;      // 新增 opticalphoton 枪
G4GenericMessenger* fMessenger = nullptr;  // 宏命令控制器
G4String fSourceMode = "gamma";            // gamma|optical|both

G4Box* fEnvelopeBox = nullptr;
const DetectorConstruction* fDetectorConstruction;
```

> 说明：`fSourceMode` 用字符串即可，简单直观，便于宏命令赋值。

---

## Step B. 修改 `src/PrimaryGeneratorAction.cc`

### B1) 在 include 区新增

```cpp
#include "G4OpticalPhoton.hh"
#include "G4GenericMessenger.hh"
```

> 你已经有 `G4SystemOfUnits.hh`，所以 `eV` 可直接用。

### B2) 在构造函数中初始化两把枪 + 命令

在你现有 `fParticleGun` 初始化后，补上：

```cpp
// 第二把枪：opticalphoton
fOpticalGun = new G4ParticleGun(1);
fOpticalGun->SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
fOpticalGun->SetParticleEnergy(2.95 * eV);  // 约 420nm

// 宏命令：/source/mode gamma|optical|both
fMessenger = new G4GenericMessenger(this, "/source/", "source control");
fMessenger->DeclareProperty("mode", fSourceMode, "gamma|optical|both");
```

### B3) 析构函数中释放

把析构函数改为：

```cpp
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fOpticalGun;
  delete fMessenger;
}
```

### B4) 修复你当前 `GeneratePrimaries()` 的作用域问题

你现在有这类写法：

```cpp
if (...) {
  G4ThreeVector Pos_Source = ...;
  G4ThreeVector Dir_Source = ...;
}
...
fParticleGun->SetParticlePosition(Pos_Source); // 这里会作用域错误
```

需要改成“先定义，再在 if 里赋值”：

```cpp
G4ThreeVector Pos_Source(0, 0, 0);
G4ThreeVector Dir_Source(0, 0, 1);

if (source_type == "Planar") {
  Pos_Source = G4ThreeVector(x0, y0, z0);
  Dir_Source = G4ThreeVector(0, 0, -z0).unit();
}
else if (source_type == "Sphere") {
  Pos_Source = r0 * URam_Pos;
  Dir_Source = URam_Dir;
}
else {
  G4cout << "Invalid source type: " << source_type << G4endl;
}
```

### B5) 在函数末尾按模式控制发射

将你原本“固定只发 `fParticleGun`”改为：

```cpp
if (fSourceMode == "gamma") {
  fParticleGun->SetParticlePosition(Pos_Source);
  fParticleGun->SetParticleMomentumDirection(Dir_Source);
  fParticleGun->GeneratePrimaryVertex(event);
}
else if (fSourceMode == "optical") {
  fOpticalGun->SetParticlePosition(Pos_Source);
  fOpticalGun->SetParticleMomentumDirection(Dir_Source);
  fOpticalGun->GeneratePrimaryVertex(event);
}
else if (fSourceMode == "both") {
  fParticleGun->SetParticlePosition(Pos_Source);
  fParticleGun->SetParticleMomentumDirection(Dir_Source);
  fParticleGun->GeneratePrimaryVertex(event);

  fOpticalGun->SetParticlePosition(Pos_Source);
  fOpticalGun->SetParticleMomentumDirection(Dir_Source);
  fOpticalGun->GeneratePrimaryVertex(event);
}
else {
  G4cerr << "Unknown /source/mode = " << fSourceMode << G4endl;
}
```

---

## 5. 宏文件怎么写（`run4mini.mac`）

在 `/run/initialize` 后加入：

```tcl
/source/mode gamma
# /source/mode optical
# /source/mode both
```

然后：

```tcl
/run/beamOn 100
```

> 解释：  
> - `/run/beamOn 100` = 跑 100 个事件  
> - 每个事件里发哪把枪 = `/source/mode` + `GeneratePrimaries()` 逻辑

---

## 6. 验证流程（建议按顺序）

1. 先设 `gamma`，确认与原结果一致（基线）
2. 再设 `optical`，看输出里是否有 `opticalphoton` track
3. 最后设 `both`，确认两种 primary 都出现

可临时加调试输出（跑通后删）：

```cpp
G4cout << "source mode = " << fSourceMode << G4endl;
```

---

## 7. 常见坑与排查

1) **`Pos_Source`/`Dir_Source` 未声明**  
- 原因：在 `if` 里定义，外面使用  
- 解法：提前定义，`if` 内赋值

2) **`/source/mode` 命令不识别**  
- 原因：没创建 `G4GenericMessenger` 或命令路径写错  
- 解法：确认构造函数里有 `new G4GenericMessenger(this, "/source/", ...)`

3) **`/gun/particle` 改了但 optical 枪不变**  
- 原因：`/gun/*` 默认只控制你原来的 gun  
- 解法：第二把枪参数要在代码里设，或后续再加 `/source/opticalEnergy` 等命令

4) **想要“物理真实闪烁”却手动打光学光子**  
- 手动打光学光子用于调试很方便  
- 真实闪烁建议靠材料+`G4OpticalPhysics` 自动产生

---

## 8. 可选增强（第二阶段）

如果你后续要更灵活，可再加这些命令：

- `/source/opticalEnergy 2.95 eV`
- `/source/opticalN 10`
- `/source/enableOptical true|false`

这样就不需要每次改 C++ 重新编译。

---

## 9. 一句话总结

你的控制链路应理解为：

`/run/beamOn X`（事件次数）  
→ `GeneratePrimaries()`（每个事件执行一次）  
→ 根据 `/source/mode` 决定调用哪把 `G4ParticleGun`

这就是“可控双源”的完整实现逻辑。


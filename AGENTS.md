# AGENTS.md
Repository guidance for autonomous coding agents in this repo.

## 1) Project Snapshot
- Stack: Geant4 C++ simulation (CMake) + Python analysis scripts.
- Main executable target: `exampleB1`.
- C++ structure: `include/*.hh`, `src/*.cc`, entrypoint `exampleB1.cc`.
- Runtime is macro-driven (`run1.mac`, `run2.mac`, `run3.mac`, `run4.mac`, `geometry.mac`, `init_vis.mac`, `vis.mac`).
- Batch helpers: `run_batch.bat`, `run_batch.sh`; visualization helper: `run_vis.bat`.
- Typical outputs go to `Results/` (with transient CSVs sometimes created in repo root first).

## 2) Cursor/Copilot Rule Files
Checked paths in this repository:
- `.cursor/rules/`: not present
- `.cursorrules`: not present
- `.github/copilot-instructions.md`: not present
If these files are added later, treat them as higher-priority instructions.

## 3) Build / Run / Test / Lint

### 3.1 Configure and Build
Windows (PowerShell/CMD):
```bash
cmake -S . -B build -DGeant4_DIR=D:/Geant4/geant4-install/lib/cmake/Geant4 -DCMAKE_PREFIX_PATH=D:/Qt2/5.15.2/msvc2019_64/lib/cmake
cmake --build build --config Release
```

POSIX:
```bash
cmake -S . -B build
cmake --build build -j
```

Ground truth from `CMakeLists.txt`:
- `find_package(Geant4 REQUIRED ui_all vis_all analysis)` is required.
- Target is `exampleB1`.
- Post-build copies macros and deploys Geant4/Qt runtime DLLs.

### 3.2 Run Commands
Windows:
```bash
build\Release\exampleB1.exe run4.mac
```

POSIX:
```bash
./build/Release/exampleB1 run4.mac
```

Convenience scripts:
```bash
run_vis.bat
run_batch.bat
./run_batch.sh
```

### 3.3 Tests and Single-Test Guidance
Current state:
- `CMakeLists.txt` has no `enable_testing()` and no `add_test(...)`.
- There is no built-in CTest suite yet.

If tests are added later:
```bash
ctest --test-dir build --output-on-failure
ctest --test-dir build -R <test_name_regex> --output-on-failure
```

Single-test equivalent today (macro smoke test):
- Run one focused macro, e.g. `exampleB1 run1.mac`.
- Use `run3.mac` for faster local checks when available.

### 3.4 Lint and Formatting
- No canonical lint/format config found (`.clang-format`, `.clang-tidy`, `ruff`, `black`, etc.).
- Do not invent a new formatting pipeline for routine tasks.
- Follow style already present in touched files.

### 3.5 Python Entrypoints (Observed)
```bash
python merge_photon_lr_perrod.py --dir <results_dir> --out merged.csv --sort
python Histo9.py --help
python Histo10_Cubic.py --help
```
Note: some older docs mention `Histo8.py`; current top-level active scripts include `Histo9.py` and `Histo10_Cubic.py`.

## 4) C++ Code Style (Observed)

### 4.1 File Organization
- Keep `.hh`/`.cc` pairs split across `include/` and `src/`.
- Keep project classes inside `namespace B1`.
- Preserve Geant4-style file header blocks where already present.

### 4.2 Includes
- In `.cc` files, include the matching local header first.
- Then include project and Geant4 headers via quotes.
- Use angle brackets for STL headers (`<filesystem>`, `<mutex>`, `<sstream>`, etc.).

### 4.3 Naming and Types
- Classes/methods: PascalCase (`DetectorConstruction`, `BeginOfRunAction`).
- Members commonly use `f` prefix (`fEdep`, `fHistoManager`, `fScoringVolume`).
- Prefer Geant4 scalar types (`G4double`, `G4int`, `G4bool`, `G4String`).
- Keep detector UI command patterns (`/detector/...`) unless task requires change.

### 4.4 Formatting and Units
- Match indentation and brace style of the touched file.
- Avoid unrelated reformatting.
- Keep explicit Geant4 unit expressions (`50 * cm`, `1.e-3 * gray`).
- Keep comments short and relevant to geometry/physics behavior.

### 4.5 Error Handling and Logging
- Use `G4cout` and `G4cerr` for runtime diagnostics.
- Avoid silent failures; use explicit early returns only when intentional.
- Prefer localized bug fixes over broad refactors in simulation-critical paths.

### 4.6 Ownership Pattern
- Existing code frequently uses explicit `new`/`delete`.
- In touched areas, follow local ownership style unless intentionally refactoring.

## 5) Python Style (Observed)
- Imports are typically stdlib first, then third-party libs.
- Prefer `pathlib.Path` for filesystem handling.
- Use `snake_case` for names and uppercase constants (`EXPECTED_COLS`).
- CSV parsing is defensive (`read_csv` fallback paths, `to_numeric(..., errors="coerce")`).
- CLI scripts use `argparse`, readable `print` diagnostics, and `sys.exit(code)`.
- Preserve local per-file style (some files typed, some legacy/minimal).

## 6) CMake and Runtime Safety Rules
- Keep executable target name `exampleB1` unless explicitly requested.
- Do not remove macro copy behavior in post-build commands.
- Preserve DLL deployment flow (`GEANT4_BIN`, optional `Qt5Core_DIR`, `CopyDlls.cmake`).
- Preserve runtime assumption: executable is run with macro files available in working directory.

## 7) Agent Working Rules
- Make focused, minimal changes tied to the request.
- Match neighboring style before applying personal/global preferences.
- Do not add dependencies unless required and justified.
- Do not commit generated artifacts (`build/`, `Results/`, transient CSV/image outputs).
- If adding tests, wire with CTest and document single-test invocation.

## 8) Validation Checklist
1. Configure/build succeeds for changed code paths.
2. Run `exampleB1` with at least one relevant macro (`run1.mac`, `run3.mac`, or `run4.mac`).
3. Confirm expected runtime outputs/logs are produced.
4. For Python changes, run modified scripts on a small sample folder.
5. Ensure generated artifacts are not accidentally staged.

## 9) Known Gaps
- `README.md` is generic/template content and not fully reliable for project specifics.
- No formal CTest suite exists in current `CMakeLists.txt`.
- No canonical lint/format pipeline is configured in-repo.

When uncertain, prioritize behavior demonstrated by `CMakeLists.txt`, macro files, and runnable scripts.

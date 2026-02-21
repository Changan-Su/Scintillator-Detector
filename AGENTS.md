# AGENTS.md

Repository guidance for autonomous coding agents.

## 1) Project Overview

- Stack: Geant4 C++ app (CMake) plus Python analysis scripts.
- Main target: `exampleB1` (`exampleB1.exe` on Windows).
- C++ layout: headers in `include/*.hh`, sources in `src/*.cc`, entrypoint `exampleB1.cc`.
- Runtime control: Geant4 macro files (`run1.mac`, `run2.mac`, `run4.mac`, `init_vis.mac`, `geometry.mac`).
- Automation scripts: `run_batch.sh`, `run_batch.bat`, `run_vis.bat`, `sub.sh`.

## 2) Cursor/Copilot Rules Check

- `.cursor/rules/`: not present.
- `.cursorrules`: not present.
- `.github/copilot-instructions.md`: not present.
- Follow repository-native conventions in this file.

## 3) Build / Run / Test / Lint

### Configure + Build

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

Notes:
- `CMakeLists.txt` requires Geant4: `find_package(Geant4 REQUIRED ui_all vis_all analysis)`.
- Post-build steps copy `.mac` files and (on Windows) Geant4/Qt DLLs to the executable directory.

### Run

Windows:

```bash
build\Release\exampleB1.exe run4.mac
```

POSIX:

```bash
./build/Release/exampleB1 run4.mac
```

Helpers:

```bash
run_vis.bat
run_batch.bat
./run_batch.sh
```

### Tests (Current State)

- No `enable_testing()` or `add_test()` currently exists in `CMakeLists.txt`.
- So there is no built-in test suite to run today.

If tests are added later, use:

```bash
ctest --test-dir build --output-on-failure
ctest --test-dir build -R <test_name_regex> --output-on-failure
```

Single-test equivalent for current repo:
- Use targeted executable + macro runs (for example `exampleB1 run1.mac`).

### Lint / Format

- No canonical lint/format configs were found (`.clang-format`, `.clang-tidy`, `ruff`, `black`, etc.).
- Do not invent lint commands; keep style consistent with surrounding files.

### Python Script Entry Points

```bash
python merge_photon_lr_perrod.py --dir <results_dir> --out merged.csv --sort
python Histo8.py <ResultsDir>
```

## 4) C++ Coding Conventions (Observed)

### File and API Structure

- Use `.hh` for headers and `.cc` for sources.
- Keep class pairs split across `include/` and `src/`.
- Keep/extend `namespace B1` patterns used in this codebase.

### Includes

- In implementation files, include the matching local header first.
- Use quoted includes for project/Geant4 headers.
- Use angle brackets for standard library headers.

### Naming

- Classes/types: PascalCase (`DetectorConstruction`, `RunAction`).
- Methods: PascalCase (`BeginOfRunAction`, `GetScoringVolume`).
- Members commonly use `f` prefix (`fScoringVolume`, `fCrystal_gap`).
- Prefer Geant4 scalar types (`G4double`, `G4int`, `G4bool`, `G4String`).

### Formatting and Units

- Match existing indentation and brace style in each file.
- Keep comments short and physics/geometry focused.
- Preserve Geant4 unit style (`50 * cm`, `56. * ns`, etc.).

### Error Handling / Logging

- Use `G4cout`/`G4cerr` for runtime diagnostics.
- Avoid silent failures; validate or clamp inputs where appropriate.
- For bugfixes, prefer minimal localized edits over refactors.

### Ownership Pattern

- Existing Geant4 lifecycle uses explicit `new`/`delete` in many places.
- Follow local ownership style in touched code to avoid mixed patterns.

## 5) Python Conventions (Observed)

- Import order: stdlib first, then third-party.
- Use `pathlib.Path` for path handling.
- Use `snake_case` for functions and variables.
- Constants are uppercase (`EXPECTED_COLS`).
- Scripts use robust CSV parsing and numeric coercion (`pd.to_numeric(..., errors="coerce")`).
- CLI scripts use `print` and `sys.exit(code)` for user-visible error handling.
- Preserve per-file style: some scripts include type hints, others do not.

## 6) CMake / Runtime Safety Rules

- Keep target name `exampleB1` stable unless explicitly requested.
- Do not remove macro copy behavior from post-build steps.
- Preserve Windows DLL deployment logic (`GEANT4_BIN` derivation and `CopyDlls.cmake` flow).

## 7) Agent Working Rules

- Make focused changes only; avoid broad cleanup.
- Match neighboring style before applying generic preferences.
- Preserve bilingual comments where already used.
- Do not add dependencies unless required by the task.
- If adding tests, wire with CTest and document single-test invocation.

## 8) Validation Checklist

1. Build/configure for the changed target succeeds.
2. Run `exampleB1` with at least one relevant macro.
3. Verify expected outputs (logs, CSVs, images, copied runtime files).
4. For Python changes, run the modified script on a small sample.
5. Avoid accidental commits of generated artifacts (`build/`, `Results/`, `Histo*_Output*`).

## 9) Known Gaps

- No canonical lint pipeline defined.
- No formal CTest suite defined today.
- `README.md` contains template text; trust source files/scripts when inconsistent.

When uncertain, prioritize behavior proven by `CMakeLists.txt`, macros, and runnable scripts.

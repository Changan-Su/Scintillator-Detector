@echo off
setlocal enabledelayedexpansion

REM ============================================================================
REM Sweep script: only source position + surface sigma
REM ============================================================================
REM Usage:
REM   1) Edit CONFIG section below
REM   2) Run: run_batch_sigma_position.bat
REM ============================================================================

REM -----------------------------
REM CONFIG
REM -----------------------------
if not defined EXE_PATH set EXE_PATH=build\Release\exampleB1.exe
if not defined GEOMETRY_MAC set GEOMETRY_MAC=geometry.mac

REM Source mode/distribution fixed for this sweep
if not defined SOURCE_MODE set SOURCE_MODE=optical
if not defined SOURCE_DISTRIBUTION set SOURCE_DISTRIBUTION=Point

REM BeamOn per run
if not defined BEAM_ON set BEAM_ON=100000

REM Optional fixed prefix before auto suffix (empty allowed)
if not defined PREFIX_BASE set PREFIX_BASE=

REM Sigma sweep
if not defined SIGMA_START set SIGMA_START=0.1
if not defined SIGMA_END set SIGMA_END=0.9
if not defined SIGMA_STEP set SIGMA_STEP=0.1

REM Position sweep (mm) for /source/fp_source x y z
if not defined X_START set X_START=10.0
if not defined X_END set X_END=10.0
if not defined X_STEP set X_STEP=1.0

if not defined Y_START set Y_START=-7.5
if not defined Y_END set Y_END=-7.5
if not defined Y_STEP set Y_STEP=1.0

if not defined Z_START set Z_START=3.0
if not defined Z_END set Z_END=3.0
if not defined Z_STEP set Z_STEP=1.0

REM Dry run: true = print only, no simulation
if not defined DRY_RUN set DRY_RUN=false

REM -----------------------------
REM ENV
REM -----------------------------
set "G4_BIN=D:\Geant4\geant4-install\bin"
set "PATH=%G4_BIN%;%PATH%"
call "D:\Geant4\geant4-install\bin\geant4.bat" >nul 2>&1

if not exist "%EXE_PATH%" (
  echo ERROR: Executable not found: %EXE_PATH%
  exit /b 1
)

if not exist "Results\" mkdir Results

set "TMP_MAC=run_tmp_sigma_position.mac"

call :generate_decimal_list SIGMA_LIST %SIGMA_START% %SIGMA_END% %SIGMA_STEP%
call :generate_decimal_list X_LIST %X_START% %X_END% %X_STEP%
call :generate_decimal_list Y_LIST %Y_START% %Y_END% %Y_STEP%
call :generate_decimal_list Z_LIST %Z_START% %Z_END% %Z_STEP%

set /a RUN_COUNT=0
echo ============================================================
echo Sweep start: sigma + source fp_source
echo ============================================================

for %%S in (%SIGMA_LIST%) do (
  for %%X in (%X_LIST%) do (
    for %%Y in (%Y_LIST%) do (
      for %%Z in (%Z_LIST%) do (
        set /a RUN_COUNT+=1

        call :sanitize_token "%%S" SIGMA_TOKEN
        call :sanitize_token "%%X" X_TOKEN
        call :sanitize_token "%%Y" Y_TOKEN
        call :sanitize_token "%%Z" Z_TOKEN

        if defined PREFIX_BASE (
          set "RUN_PREFIX=!PREFIX_BASE!_S!SIGMA_TOKEN!_X!X_TOKEN!_Y!Y_TOKEN!_Z!Z_TOKEN!"
        ) else (
          set "RUN_PREFIX=S!SIGMA_TOKEN!_X!X_TOKEN!_Y!Y_TOKEN!_Z!Z_TOKEN!"
        )

        echo [!RUN_COUNT!] sigma=%%S  pos=^(%%X, %%Y, %%Z^)  prefix=!RUN_PREFIX!

        (
          echo /control/execute %GEOMETRY_MAC%
          echo /control/verbose 0
          echo /run/verbose 1
          echo /event/verbose 0
          echo /tracking/verbose 0
          echo /run/printProgress 100
          echo /vis/disable
          echo /detector/surfaceSigma %%S
          echo /results/prefix !RUN_PREFIX!
          echo /run/initialize
          echo /source/mode %SOURCE_MODE%
          echo /source/distribution %SOURCE_DISTRIBUTION%
          echo /source/fp_source %%X %%Y %%Z
          echo /run/beamOn %BEAM_ON%
        ) > "%TMP_MAC%"

        if /I "%DRY_RUN%"=="true" (
          echo     DRY_RUN=true, skip simulation.
        ) else (
          "%EXE_PATH%" "%TMP_MAC%"
          if errorlevel 1 (
            echo     ERROR: simulation failed.
          ) else (
            for /f "delims=" %%D in ('dir /b /ad /o-d "Results"') do (
              if not defined LAST_DIR set "LAST_DIR=%%D"
            )
            if defined LAST_DIR (
              echo     Latest: Results\!LAST_DIR!
              set "LAST_DIR="
            )
          )
        )
        echo.
      )
    )
  )
)

del "%TMP_MAC%" >nul 2>&1
echo ============================================================
echo Sweep done. Total runs: %RUN_COUNT%
echo ============================================================
exit /b 0

:sanitize_token
setlocal
set "V=%~1"
set "V=%V:.=p%"
set "V=%V:-=m%"
endlocal & set "%~2=%V%"
exit /b 0

:generate_decimal_list
setlocal enabledelayedexpansion
set "LIST="
for /f "delims=" %%V in ('powershell -NoProfile -Command "$s=[double]('%~2');$e=[double]('%~3');$st=[double]('%~4');if($st -eq 0){exit 1};if($st -gt 0){for($v=$s;$v -le $e+1e-9;$v+=$st){'{0:0.0}' -f $v}}else{for($v=$s;$v -ge $e-1e-9;$v+=$st){'{0:0.0}' -f $v}}"') do (
  set "LIST=!LIST! %%V"
)
endlocal & set "%~1=%LIST%"
exit /b 0

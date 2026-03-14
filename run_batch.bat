@echo off
REM ============================================================================
REM Parameter Sweep Batch Script for Geant4 Scintillator Detector
REM Author: Auto-generated
REM Date: 2026-01-31
REM ============================================================================
REM
REM This script automatically loops through geometry parameters, runs
REM simulations, and stores results in separate folders.
REM
REM Usage: Edit configuration section below, then run: run_batch.bat
REM ============================================================================

setlocal enabledelayedexpansion

REM ============================================================================
REM CONFIGURATION SECTION - Edit these values
REM ============================================================================

REM --- Parameter Loop Configuration (true/false) ---
set LOOP_ARRAY_NX=false
set LOOP_ARRAY_NY=false
set LOOP_ARRAY_NZ=false
set LOOP_CRYSTAL_GAP=false
set LOOP_CRYSTAL_SIZE=false
set LOOP_CRYSTAL_SIZE_Y=false
set LOOP_FILLTER_RATIO_Y=false
set LOOP_FILLTER_RATIO_Z=false
set LOOP_FILLTER_POS_RATIO_Y=false
set LOOP_FILLTER_POS_RATIO_Z=false
set LOOP_SURFACE_SIGMA=false

REM --- Parameter Range Configuration ---
REM Array Nx (crystal count in x direction)
set NX_START=11
set NX_END=11
set NX_STEP=1

REM Array Ny (crystal count in y direction)
set NY_START=1
set NY_END=3
set NY_STEP=1

REM Array Nz (crystal count in z direction)
set NZ_START=7
set NZ_END=7
set NZ_STEP=1

REM Crystal Gap (in mm, 0 means continuous no-gap packing)
set GAP_START=0.0
set GAP_END=0.0
set GAP_STEP=0.1

REM Crystal Size (x/z dimension, in mm)
set SIZE_START=3
set SIZE_END=3
set SIZE_STEP=1

REM Crystal Size Y (y dimension, in mm)
set SIZE_Y_START=3
set SIZE_Y_END=3
set SIZE_Y_STEP=1

REM Fillter Ratio Y (0-1, size ratio in Y)
set FILLTER_RATIO_Y_START=0.3
set FILLTER_RATIO_Y_END=0.3
set FILLTER_RATIO_Y_STEP=0.1

REM Fillter Ratio Z (0-1, size ratio in Z)
set FILLTER_RATIO_Z_START=1.0
set FILLTER_RATIO_Z_END=1.0
set FILLTER_RATIO_Z_STEP=0.1

REM Fillter Position Ratio Y (0-1, position ratio in Y)
set FILLTER_POS_RATIO_Y_START=0.0
set FILLTER_POS_RATIO_Y_END=0.9
set FILLTER_POS_RATIO_Y_STEP=0.1

REM Fillter Position Ratio Z (0-1, position ratio in Z)
set FILLTER_POS_RATIO_Z_START=0.0
set FILLTER_POS_RATIO_Z_END=0.0
set FILLTER_POS_RATIO_Z_STEP=0.1

REM Surface Sigma (0-1, crystal optical surface roughness)
set SIGMA_START=0.5
set SIGMA_END=0.7
set SIGMA_STEP=0.1

REM --- Default Values (used when parameter is not looping) ---
set DEFAULT_NX=9
set DEFAULT_NY=1
set DEFAULT_NZ=1
set DEFAULT_GAP=0.0
set DEFAULT_SIZE=3
set DEFAULT_SIZE_Y=3
set DEFAULT_FILLTER_RATIO_Y=0.3
set DEFAULT_FILLTER_RATIO_Z=1.0
set DEFAULT_FILLTER_POS_RATIO_Y=0.7
set DEFAULT_FILLTER_POS_RATIO_Z=0.0
set DEFAULT_SURFACE_SIGMA=0.5

REM --- Run Configuration ---
set RUN_MACRO=run4.mac
set EXE_PATH=build\Release\exampleB1.exe

REM --- Folder Naming Configuration (true/false) ---
set NAME_INCLUDE_NX=false
set NAME_INCLUDE_NY=false
set NAME_INCLUDE_NZ=false
set NAME_INCLUDE_GAP=false
set NAME_INCLUDE_SIZE=false
set NAME_INCLUDE_SIZE_Y=false
set NAME_INCLUDE_FILLTER_RATIO_Y=false
set NAME_INCLUDE_FILLTER_RATIO_Z=false
set NAME_INCLUDE_FILLTER_POS_RATIO_Y=true
set NAME_INCLUDE_FILLTER_POS_RATIO_Z=false
set NAME_INCLUDE_SURFACE_SIGMA=false

REM ============================================================================
REM END OF CONFIGURATION - Do not edit below unless you know what you're doing
REM ============================================================================

echo ============================================================================
echo Parameter Sweep Batch Script
echo ============================================================================
echo.

REM Set up Geant4 environment
set "G4_BIN=D:\Geant4\geant4-install\bin"
set "PATH=%G4_BIN%;%PATH%"
call "D:\Geant4\geant4-install\bin\geant4.bat" >nul 2>&1

REM Check if executable exists
if not exist "%EXE_PATH%" (
    echo ERROR: Executable not found at %EXE_PATH%
    echo Please build the project first.
    pause
    exit /b 1
)

REM Create Results directory
if not exist "Results\" mkdir Results

REM Initialize loop counter
set /a LOOP_COUNT=0

REM Start parameter loops
echo Starting parameter sweep...
echo.

REM Loop through Nx
set "NX_LIST="
if "%LOOP_ARRAY_NX%"=="true" (
    for /L %%N in (%NX_START%,%NX_STEP%,%NX_END%) do (
        set "NX_LIST=!NX_LIST! %%N"
    )
) else (
    set "NX_LIST=%DEFAULT_NX%"
)

REM Loop through Ny
set "NY_LIST="
if "%LOOP_ARRAY_NY%"=="true" (
    for /L %%N in (%NY_START%,%NY_STEP%,%NY_END%) do (
        set "NY_LIST=!NY_LIST! %%N"
    )
) else (
    set "NY_LIST=%DEFAULT_NY%"
)

REM Loop through Nz
set "NZ_LIST="
if "%LOOP_ARRAY_NZ%"=="true" (
    for /L %%N in (%NZ_START%,%NZ_STEP%,%NZ_END%) do (
        set "NZ_LIST=!NZ_LIST! %%N"
    )
) else (
    set "NZ_LIST=%DEFAULT_NZ%"
)

REM Loop through Gap (using decimal values)
set "GAP_LIST="
if "%LOOP_CRYSTAL_GAP%"=="true" (
    call :generate_decimal_list GAP_LIST %GAP_START% %GAP_END% %GAP_STEP%
) else (
    set "GAP_LIST=%DEFAULT_GAP%"
)

REM Loop through Size
set "SIZE_LIST="
if "%LOOP_CRYSTAL_SIZE%"=="true" (
    call :generate_decimal_list SIZE_LIST %SIZE_START% %SIZE_END% %SIZE_STEP%
) else (
    set "SIZE_LIST=%DEFAULT_SIZE%"
)

REM Loop through Size Y
set "SIZE_Y_LIST="
if "%LOOP_CRYSTAL_SIZE_Y%"=="true" (
    call :generate_decimal_list SIZE_Y_LIST %SIZE_Y_START% %SIZE_Y_END% %SIZE_Y_STEP%
) else (
    set "SIZE_Y_LIST=%DEFAULT_SIZE_Y%"
)

REM Loop through Fillter Ratio Y
set "FILLTER_RATIO_Y_LIST="
if "%LOOP_FILLTER_RATIO_Y%"=="true" (
    call :generate_decimal_list FILLTER_RATIO_Y_LIST %FILLTER_RATIO_Y_START% %FILLTER_RATIO_Y_END% %FILLTER_RATIO_Y_STEP%
) else (
    set "FILLTER_RATIO_Y_LIST=%DEFAULT_FILLTER_RATIO_Y%"
)

REM Loop through Fillter Ratio Z
set "FILLTER_RATIO_Z_LIST="
if "%LOOP_FILLTER_RATIO_Z%"=="true" (
    call :generate_decimal_list FILLTER_RATIO_Z_LIST %FILLTER_RATIO_Z_START% %FILLTER_RATIO_Z_END% %FILLTER_RATIO_Z_STEP%
) else (
    set "FILLTER_RATIO_Z_LIST=%DEFAULT_FILLTER_RATIO_Z%"
)

REM Loop through Fillter Position Ratio Y
set "FILLTER_POS_RATIO_Y_LIST="
if "%LOOP_FILLTER_POS_RATIO_Y%"=="true" (
    call :generate_decimal_list FILLTER_POS_RATIO_Y_LIST %FILLTER_POS_RATIO_Y_START% %FILLTER_POS_RATIO_Y_END% %FILLTER_POS_RATIO_Y_STEP%
) else (
    set "FILLTER_POS_RATIO_Y_LIST=%DEFAULT_FILLTER_POS_RATIO_Y%"
)

REM Loop through Fillter Position Ratio Z
set "FILLTER_POS_RATIO_Z_LIST="
if "%LOOP_FILLTER_POS_RATIO_Z%"=="true" (
    call :generate_decimal_list FILLTER_POS_RATIO_Z_LIST %FILLTER_POS_RATIO_Z_START% %FILLTER_POS_RATIO_Z_END% %FILLTER_POS_RATIO_Z_STEP%
) else (
    set "FILLTER_POS_RATIO_Z_LIST=%DEFAULT_FILLTER_POS_RATIO_Z%"
)

REM Loop through Surface Sigma
set "SIGMA_LIST="
if "%LOOP_SURFACE_SIGMA%"=="true" (
    call :generate_decimal_list SIGMA_LIST %SIGMA_START% %SIGMA_END% %SIGMA_STEP%
) else (
    set "SIGMA_LIST=%DEFAULT_SURFACE_SIGMA%"
)

REM Nested loops through all parameter combinations
for %%5 in (%SIGMA_LIST%) do (
    for %%x in (%NX_LIST%) do (
        for %%y in (%NY_LIST%) do (
        for %%z in (%NZ_LIST%) do (
            for %%g in (%GAP_LIST%) do (
                for %%s in (%SIZE_LIST%) do (
                    for %%t in (%SIZE_Y_LIST%) do (
                        for %%1 in (%FILLTER_RATIO_Y_LIST%) do (
                            for %%2 in (%FILLTER_RATIO_Z_LIST%) do (
                                for %%3 in (%FILLTER_POS_RATIO_Y_LIST%) do (
                                    for %%4 in (%FILLTER_POS_RATIO_Z_LIST%) do (
                                        set /a LOOP_COUNT+=1
                                        
                                        REM Build folder name
                                        set "FOLDER_NAME=!LOOP_COUNT!"
                                        if !LOOP_COUNT! LSS 10 set "FOLDER_NAME=00!LOOP_COUNT!"
                                        if !LOOP_COUNT! GEQ 10 if !LOOP_COUNT! LSS 100 set "FOLDER_NAME=0!LOOP_COUNT!"
                                        
                                        if "%NAME_INCLUDE_NX%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_Nx%%x"
                                        if "%NAME_INCLUDE_NY%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_Ny%%y"
                                        if "%NAME_INCLUDE_NZ%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_Nz%%z"
                                        if "%NAME_INCLUDE_GAP%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_Gap%%g"
                                        if "%NAME_INCLUDE_SIZE%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_Size%%s"
                                        if "%NAME_INCLUDE_SIZE_Y%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_SizeY%%t"
                                        if "%NAME_INCLUDE_FILLTER_RATIO_Y%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_FRY%%1"
                                        if "%NAME_INCLUDE_FILLTER_RATIO_Z%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_FRZ%%2"
                                        if "%NAME_INCLUDE_FILLTER_POS_RATIO_Y%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_FPRY%%3"
                                        if "%NAME_INCLUDE_FILLTER_POS_RATIO_Z%"=="true" set "FOLDER_NAME=!FOLDER_NAME!_FPRZ%%4"
                                        if "%NAME_INCLUDE_SURFACE_SIGMA%"=="true" (
                                            set "SIGMA_STR=%%5"
                                            set "SIGMA_STR=!SIGMA_STR:.=p!"
                                            set "FOLDER_NAME=!FOLDER_NAME!_Sigma!SIGMA_STR!"
                                        )
                                        
                                        echo [!LOOP_COUNT!] Running: Nx=%%x Ny=%%y Nz=%%z Gap=%%g Size=%%s SizeY=%%t FillterY=%%1 FillterZ=%%2 FPosY=%%3 FPosZ=%%4 Sigma=%%5
                                        echo    Output: Results\^<timestamp+params^>
                                        
                                        REM Generate geometry.mac
                                        call :generate_geometry_mac %%x %%y %%z %%g %%s %%t %%1 %%2 %%3 %%4 %%5
                                        
                                        REM Run simulation
                                        "%EXE_PATH%" "%RUN_MACRO%"
                                        if errorlevel 1 (
                                            echo    ERROR: Simulation failed for this parameter set.
                                        ) else (
                                            set "LATEST_RESULT_DIR="
                                            for /f "delims=" %%D in ('dir /b /ad /o-d "Results"') do (
                                                if not defined LATEST_RESULT_DIR set "LATEST_RESULT_DIR=%%D"
                                            )
                                            if defined LATEST_RESULT_DIR (
                                                echo    Latest output dir: Results\!LATEST_RESULT_DIR!
                                            )
                                        )
                                        
                                        echo    Completed.
                                        echo.
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
)

echo ============================================================================
echo Parameter sweep completed!
echo Total runs: %LOOP_COUNT%
echo Results saved in: Results\
echo ============================================================================
pause
exit /b 0

REM ============================================================================
REM Subroutines
REM ============================================================================

:generate_geometry_mac
REM Generate geometry.mac with parameters: Nx Ny Nz Gap Size SizeY FillterRatioY FillterRatioZ FillterPosRatioY FillterPosRatioZ SurfaceSigma
REM Workaround for %10+ parameters: store them first
set "GM_NX=%1"
set "GM_NY=%2"
set "GM_NZ=%3"
set "GM_GAP=%4"
set "GM_SIZE=%5"
set "GM_SIZEY=%6"
set "GM_FRY=%7"
set "GM_FRZ=%8"
set "GM_FPRY=%9"
shift
set "GM_FPRZ=%9"
shift
set "GM_SIGMA=%9"

(
echo # Geometry macro: crystal array and gap ^(run before /run/initialize^)
echo # Auto-generated by run_batch.bat
echo #
echo # Crystal array ^(nx, ny, nz^)
echo /detector/arrayNx %GM_NX%
echo /detector/arrayNy %GM_NY%
echo /detector/arrayNz %GM_NZ%
echo #
echo # Crystal gap in mm ^(0 means continuous no-gap packing^)
echo /detector/crystalGap %GM_GAP% mm
echo #
echo # Single crystal size ^(x/z^), in mm
echo /detector/crystalSize %GM_SIZE% mm
echo #
echo # Single crystal size in y, in mm
echo /detector/crystalSizeY %GM_SIZEY% mm
echo #
echo # Fillter ^(gap filler^) parameters ^(ratios 0-1^)
echo /detector/fillterRatioY %GM_FRY%
echo /detector/fillterRatioZ %GM_FRZ%
echo /detector/fillterPosRatioY %GM_FPRY%
echo /detector/fillterPosRatioZ %GM_FPRZ%
echo #
echo # Crystal optical surface roughness ^(0-1^)
echo /detector/surfaceSigma %GM_SIGMA%
) > geometry.mac
exit /b 0

:generate_decimal_list
REM Generate list of decimal values: VAR_NAME START END STEP
setlocal enabledelayedexpansion
set "LIST="
set "START_INT=%~2"
set "END_INT=%~3"
set "STEP_INT=%~4"

REM Convert to integer arithmetic (multiply by 10 to handle one decimal place)
REM Remove decimal point and strip leading zeros to avoid octal interpretation
set "START_STR=!START_INT:.=!"
set "END_STR=!END_INT:.=!"
set "STEP_STR=!STEP_INT:.=!"

REM Strip leading zeros by forcing decimal base (add 10# prefix or use string manipulation)
for /f "tokens=* delims=0" %%a in ("!START_STR!") do set "START_STR=%%a"
for /f "tokens=* delims=0" %%a in ("!END_STR!") do set "END_STR=%%a"
for /f "tokens=* delims=0" %%a in ("!STEP_STR!") do set "STEP_STR=%%a"

REM Handle edge case where value is exactly 0 or 00
if "!START_STR!"=="" set "START_STR=0"
if "!END_STR!"=="" set "END_STR=0"
if "!STEP_STR!"=="" set "STEP_STR=0"

set /a "START_INT_X10=!START_STR!"
set /a "END_INT_X10=!END_STR!"
set /a "STEP_INT_X10=!STEP_STR!"

REM Handle case where input doesn't have decimal point
if "!START_INT!"=="!START_INT:.=!" set /a "START_INT_X10=!START_INT! * 10"
if "!END_INT!"=="!END_INT:.=!" set /a "END_INT_X10=!END_INT! * 10"
if "!STEP_INT!"=="!STEP_INT:.=!" set /a "STEP_INT_X10=!STEP_INT! * 10"

for /L %%N in (!START_INT_X10!,!STEP_INT_X10!,!END_INT_X10!) do (
    set /a "WHOLE=%%N / 10"
    set /a "FRAC=%%N %% 10"
    set "LIST=!LIST! !WHOLE!.!FRAC!"
)

endlocal & set "%~1=%LIST%"
exit /b 0

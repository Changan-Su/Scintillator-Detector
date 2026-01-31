@echo off
REM Run exampleB1 with Qt visualization
REM Execute from project root: run_vis.bat

set "G4_BIN=D:\Geant4\geant4-install\bin"
set "PATH=%G4_BIN%;%PATH%"

set "QT_BIN=D:\Qt2\5.15.2\msvc2019_64\bin"
set "PATH=%QT_BIN%;%PATH%"

call "D:\Geant4\geant4-install\bin\geant4.bat"

cd /d "%~dp0"
REM Run exe from project root so geometry.mac and init_vis.mac are read from here (no rebuild needed after editing .mac)
if exist "build\Release\exampleB1.exe" (
  build\Release\exampleB1.exe
) else if exist "D:\Geant4\geant4-install\share\Geant4\examples\basic\B1\Release\exampleB1.exe" (
  "D:\Geant4\geant4-install\share\Geant4\examples\basic\B1\Release\exampleB1.exe"
) else (
  echo Build folder not found. Run: mkdir build ^&^& cd build ^&^& cmake .. ^&^& cmake --build . --config Release
)

pause

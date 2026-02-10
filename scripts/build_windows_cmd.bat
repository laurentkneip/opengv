@echo off
REM Build script for Windows using MSVC (cmd.exe)
REM Usage: scripts\build_windows_cmd.bat

setlocal enabledelayedexpansion

cd /d "%~dp0\.."

echo ========================================
echo OpenGV Windows Build (MSVC)
echo ========================================

REM Activate conda environment
call conda activate opengv
if %errorlevel% neq 0 (
    echo ERROR: Failed to activate conda environment
    exit /b 1
)
echo [OK] Conda environment activated

echo.
echo Step 1: Installing Python dependencies...
pip install numpy scipy pytest >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Failed to install Python dependencies
    exit /b 1
)
echo [OK] Python dependencies installed

echo.
echo Step 2: Configuring CMake...
if exist build rmdir /s /q build
mkdir build
cd build
cmake .. -DBUILD_PYTHON=ON -DBUILD_TESTS=ON -A x64 -DCMAKE_CXX_FLAGS="/MP /EHsc"
if %errorlevel% neq 0 (
    echo ERROR: CMake configuration failed
    exit /b 1
)
echo [OK] CMake configured successfully

echo.
echo Step 3: Building (Release, parallel)...
cmake --build . --config Release --parallel
if %errorlevel% neq 0 (
    echo ERROR: Build failed
    exit /b 1
)
echo [OK] Build completed successfully

echo.
echo Step 4: Running C++ tests (verbose for debugging)...
ctest -C Release --output-on-failure --verbose
if %errorlevel% neq 0 (
    echo [WARN] Some C++ tests failed
) else (
    echo [OK] All C++ tests passed
)

echo.
echo Step 5: Checking outputs...
if exist lib\Release\opengv.lib (
    echo [OK] Found: lib\Release\opengv.lib
) else (
    echo [WARN] opengv.lib not found in lib\Release
)

if exist lib\Release\pyopengv.pyd (
    echo [OK] Found: lib\Release\pyopengv.pyd
) else (
    echo [WARN] pyopengv.pyd not found in lib\Release
)

echo.
echo Step 6: Testing Python module...
cd ..
set PYTHONPATH=%CD%\build\lib\Release
pytest python\ -v
if %errorlevel% neq 0 (
    echo [WARN] Some Python tests failed
) else (
    echo [OK] All Python tests passed
)

echo.
echo ========================================
echo Windows build completed!
echo ========================================

endlocal

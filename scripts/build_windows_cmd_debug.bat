@echo off
REM Debug build script for Windows using MSVC (cmd.exe)
REM This script builds in Debug mode with all sanitizers and verbose output
REM Usage: scripts\build_windows_cmd_debug.bat

setlocal enabledelayedexpansion

cd /d "%~dp0\.."

echo ========================================
echo OpenGV Windows Debug Build (MSVC)
echo ========================================
echo.
echo This build includes:
echo  - Debug symbols
echo  - No optimizations
echo  - Address sanitizer
echo  - Verbose test output
echo  - Stack traces on crashes
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
pip install --no-build-isolation -e ".[dev]" >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Failed to install Python dependencies
    exit /b 1
)
echo [OK] Python dependencies installed

echo.
echo Step 2: Configuring CMake (Debug mode with sanitizers)...
if exist build_debug rmdir /s /q build_debug
mkdir build_debug
cd build_debug

REM Configure with Debug mode and address sanitizer
cmake .. ^
    -DBUILD_PYTHON=ON ^
    -DBUILD_TESTS=ON ^
    -A x64 ^
    -DCMAKE_BUILD_TYPE=Debug ^
    -DCMAKE_CXX_FLAGS="/Od /Zi /DEBUG:FULL /fsanitize=address /MDd"

if %errorlevel% neq 0 (
    echo ERROR: CMake configuration failed
    exit /b 1
)
echo [OK] CMake configured successfully

echo.
echo Step 3: Building (Debug mode, single-threaded for better error messages)...
cmake --build . --config Debug -j 1
if %errorlevel% neq 0 (
    echo ERROR: Build failed
    exit /b 1
)
echo [OK] Build completed successfully

echo.
echo Step 4: Checking outputs...
if exist lib\Debug\opengv.lib (
    echo [OK] Found: lib\Debug\opengv.lib
) else (
    echo [WARN] opengv.lib not found in lib\Debug
)

if exist lib\Debug\pyopengv.pyd (
    echo [OK] Found: lib\Debug\pyopengv.pyd
) else (
    echo [WARN] pyopengv.pyd not found in lib\Debug
)

echo.
echo Step 5: Running C++ tests (VERBOSE with Debug symbols)...
echo ============================================================
echo Test output will show exact line where crashes occur
echo ============================================================
echo.

REM Run with maximum verbosity and stop on first failure
ctest -C Debug ^
    --output-on-failure ^
    --verbose ^
    --debug ^
    --stop-on-failure

if %errorlevel% neq 0 (
    echo.
    echo ============================================================
    echo [ERROR] Test failed with segfault or assertion
    echo ============================================================
    echo.
    echo The crash occurred in Debug mode with:
    echo  - Full debug symbols loaded
    echo  - Address sanitizer active
    echo  - Stack trace available above
    echo.
    echo To debug further:
    echo  1. Look at the last printed line before crash
    echo  2. Check for ASAN error messages above
    echo  3. Run individual test with: build_debug\Debug\test_name.exe
    echo  4. Use Visual Studio debugger: devenv build_debug\opengv.sln
    echo.
    pause
    exit /b 1
) else (
    echo [OK] All C++ tests passed
)

echo.
echo Step 6: Testing Python module...
cd ..
set PYTHONPATH=%CD%\build_debug\lib\Debug
pytest python\ -v
if %errorlevel% neq 0 (
    echo [WARN] Some Python tests failed
) else (
    echo [OK] All Python tests passed
)

echo.
echo ========================================
echo Windows Debug build completed!
echo ========================================
echo.
echo Output location: build_debug\
echo.
echo To debug a specific test manually:
echo   cd build_debug\Debug
echo   test_relative_pose.exe
echo.
echo To debug in Visual Studio:
echo   devenv build_debug\opengv.sln
echo   Set test as startup project
echo   Press F5 to debug
echo.

endlocal

@echo off
REM Test script to validate conda environment setup for OpenGV on Windows
REM Usage: scripts\test_conda_env.bat

setlocal enabledelayedexpansion

REM Change to project root directory
cd /d "%~dp0\.."

echo Testing OpenGV Conda Environment
echo =================================

REM Check if conda is available
where conda >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: conda is not installed or not in PATH
    exit /b 1
)

REM Check if environment exists
conda env list | findstr /C:"opengv " >nul
if %errorlevel% neq 0 (
    echo ERROR: 'opengv' conda environment does not exist
    echo Please create it first: conda env create -f environment.yml
    exit /b 1
)

echo [OK] Conda environment 'opengv' exists

REM Activate environment
call conda activate opengv
if %errorlevel% neq 0 (
    echo ERROR: Failed to activate conda environment
    exit /b 1
)

echo [OK] Activated conda environment

REM Check Python
echo Checking Python...
python --version
if %errorlevel% neq 0 (
    echo ERROR: Python not found
    exit /b 1
)
python -c "import sys; assert sys.version_info >= (3, 7)" >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Python 3.7+ required
    exit /b 1
)
echo [OK] Python version OK

REM Check required packages
echo Checking required packages...
set MISSING=0

python -c "import numpy" >nul 2>&1
if %errorlevel% equ 0 (
    echo [OK] numpy
) else (
    echo [FAIL] numpy MISSING
    set MISSING=1
)

python -c "import pybind11" >nul 2>&1
if %errorlevel% equ 0 (
    echo [OK] pybind11
) else (
    echo [FAIL] pybind11 MISSING
    set MISSING=1
)

if !MISSING! equ 1 (
    echo ERROR: Missing required packages
    exit /b 1
)

REM Check CMake
echo Checking CMake...
where cmake >nul 2>&1
if %errorlevel% equ 0 (
    cmake --version
    echo [OK] CMake found
) else (
    echo ERROR: CMake not found
    exit /b 1
)

REM Check Eigen (via conda or system)
echo Checking Eigen...
conda list eigen | findstr /C:"eigen" >nul
if %errorlevel% equ 0 (
    echo [OK] Eigen found via conda
) else if exist "%CONDA_PREFIX%\include\eigen3" (
    echo [OK] Eigen headers found
) else if exist "%CONDA_PREFIX%\include\eigen" (
    echo [OK] Eigen headers found
) else (
    echo [WARN] Eigen headers not found in expected location (may still work if in system path)
)

REM Test basic build configuration
echo.
echo Testing CMake configuration...
if exist build_test rmdir /s /q build_test
mkdir build_test
cd build_test

cmake .. -DBUILD_PYTHON=ON -DBUILD_TESTS=OFF -A x64 >cmake_output.log 2>&1
if %errorlevel% equ 0 (
    echo [OK] CMake configuration successful
    
    REM Check if Python bindings target was created
    findstr /C:"pyopengv" cmake_output.log >nul
    if %errorlevel% equ 0 (
        echo [OK] Python bindings target configured
    )
    
    cd ..
    rmdir /s /q build_test
) else (
    echo ERROR: CMake configuration failed
    type cmake_output.log
    cd ..
    rmdir /s /q build_test
    exit /b 1
)

echo.
echo =================================
echo All tests passed! [OK]
echo The conda environment is ready for building OpenGV.
echo.
echo To build OpenGV:
echo   1. mkdir build ^&^& cd build
echo   2. cmake .. -DBUILD_PYTHON=ON -A x64
echo   3. cmake --build . --config Release

endlocal

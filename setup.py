#!/usr/bin/env python3
"""
Setup script for pyopengv - Python bindings for OpenGV library.

This setup.py provides options to:
1. Build and install locally as a wheel
2. Upload to PyPI

Usage:
    # Install locally (development mode)
    pip install -e .
    
    # Build wheel locally
    pip install build
    python -m build
    
    # Install from wheel
    pip install dist/pyopengv-*.whl
    
    # Upload to PyPI (requires twine)
    pip install twine
    python -m build
    twine upload dist/*
    
    # Upload to TestPyPI first (recommended)
    twine upload --repository testpypi dist/*
"""

import os
import sys
import subprocess
import platform
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """Extension that is built using CMake."""
    
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """Custom build extension that runs CMake."""
    
    def build_extension(self, ext):
        if not isinstance(ext, CMakeExtension):
            super().build_extension(ext)
            return
            
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        # Required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        # CMake configuration arguments
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            "-DBUILD_TESTS=OFF",
            "-DBUILD_PYTHON=ON",
        ]

        # Build type
        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        # Platform-specific configuration
        if platform.system() == "Windows":
            cmake_args += [
                f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}",
            ]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += [f"-DCMAKE_BUILD_TYPE={cfg}"]
            # Parallel build
            import multiprocessing
            num_jobs = multiprocessing.cpu_count()
            build_args += ["--", f"-j{num_jobs}"]

        env = os.environ.copy()
        env["CXXFLAGS"] = f'{env.get("CXXFLAGS", "")} -DVERSION_INFO=\\"{self.distribution.get_version()}\\"'

        # Create build directory
        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        # Run CMake configure
        print(f"Running CMake configuration in {build_temp}")
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args,
            cwd=build_temp,
            env=env
        )

        # Run CMake build
        print(f"Building extension with CMake")
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args,
            cwd=build_temp
        )


def main():
    """Main setup function."""
    
    # Read version from file if available
    version = "1.0.0"

    # Read long description from README
    readme_file = Path(__file__).parent / "README.md"
    long_description = ""
    if readme_file.exists():
        long_description = readme_file.read_text(encoding="utf-8")
        long_description_content_type = "text/markdown"
    else:
        long_description_content_type = "text/plain"
    
    setup(
        name="pyopengv",
        version=version,
        author="Laurent Kneip",
        author_email="kneip.laurent@gmail.com",
        description="Python bindings for OpenGV - geometric vision library",
        long_description=long_description,
        long_description_content_type=long_description_content_type,
        ext_modules=[CMakeExtension("pyopengv")],
        cmdclass={"build_ext": CMakeBuild},
        zip_safe=False,
        python_requires=">=3.7",
    )


if __name__ == "__main__":
    main()

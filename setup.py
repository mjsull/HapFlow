# cx_Freeze setup file

import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["pysam"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "HapFlow",
        version = "1.1.2",
        description = "Visualising haplotypes in sequencing data.",
        options = {"build_exe": build_exe_options},
        executables = [Executable("HapFlow.py", base=base)])
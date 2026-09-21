# ----------------------------------------------------------------------
# Config file for LLVM (clang/flang)
# Note: flang only accepts -std=f2018 (no f2008), so no -std flag is set.
# ----------------------------------------------------------------------

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpifort
CC   = mpicc

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS = -fPIC -fdefault-real-8 -O2
FF90_FLAGS = ${FF77_FLAGS}
C_FLAGS    = -fPIC -O2

# ------- Define Archiver and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS = -fPIC

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config
F2PY = f2py

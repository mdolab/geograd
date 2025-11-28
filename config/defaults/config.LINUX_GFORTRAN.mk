# ----------------------------------------------------------------------
# Config file for gfortran
# ----------------------------------------------------------------------

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpifort
CC   = mpicc

# ------- Define complexify inlcude and linker flags -------------------------
#COMPLEXIFY_INCLUDE_FLAGS=-I$(COMPLEXIFY_DIR)/include
#COMPLEXIFY_LINKER_FLAGS=-L$(COMPLEXIFY_DIR)/lib -lcomplexify

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS = -fPIC -fdefault-real-8 -O2
FF90_FLAGS = ${FF77_FLAGS} -std=f2008
C_FLAGS    = -fPIC -O2

# ------- Define Archiver and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS = -fPIC


# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py

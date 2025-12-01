# ----------------------------------------------------------------------
# Config file for Intel ifort
# ----------------------------------------------------------------------

# ------- Define the MPI Compilers--------------------------------------
ifdef I_MPI_ROOT # Using Intel MPI
  # Note that ";" is there to avoid make shell optimization, otherwise the shell command may fail
  ICC_EXISTS := $(shell command -v icc;)

  ifdef ICC_EXISTS
    # icc only exists on older Intel versions
    # Assume that we want to use the old compilers
    FF90 = mpiifort
    CC = mpiicc
  else
    # Use the new compilers
    FF90 = mpiifx
    CC = mpiicx
  endif
else # Using HPE MPI
  FF90 = ifort -lmpi
  CC   = icc -lmpi
endif

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS = -fPIC -r8 -O2
FF90_FLAGS = ${FF77_FLAGS} #-stand f08
C_FLAGS    = -fPIC -O2

# ------- Define Archiver and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS = -fPIC -nofor-main

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config
F2PY = f2py

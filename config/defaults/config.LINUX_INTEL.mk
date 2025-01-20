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
    FCOMPILER = mpiifort
  else
    # Use the new compilers
    FCOMPILER = mpiifx
  endif
else # Using HPE MPI
  FCOMPILER = ifort -lmpi
endif

FCOMPILER_ALL_FLAGS=$(FCOMPILER) -O3
F2PY_INCLUDES=-I$(MPI_INSTALL_DIR)/include
F2PY_ALL_FLAGS=--fcompiler=intelem --f90exec=$(FCOMPILER) --f77exec=$(FCOMPILER) $(F2PY_INCLUDES) --opt='-O3'
F2PY=python3 -m numpy.f2py
                                  
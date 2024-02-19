# ----------------------------------------------------------------------
# Config file for Intel ifort
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
ifdef I_MPI_ROOT # Using Intel MPI
  FCOMPILER = mpiifort
else # Using HPE MPI
  FCOMPILER = ifort -lmpi
endif

FCOMPILER_ALL_FLAGS=$(FCOMPILER) -O3
F2PY_INCLUDES=-I$(MPI_INSTALL_DIR)/include
F2PY_ALL_FLAGS=--fcompiler=intelem --f90exec=$(FCOMPILER) --f77exec=$(FCOMPILER) $(F2PY_INCLUDES) --opt='-O3'
F2PY=python3 -m numpy.f2py
                                  
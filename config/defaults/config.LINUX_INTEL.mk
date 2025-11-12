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
  else
    # Use the new compilers
    FF90 = mpiifx
  endif
else # Using HPE MPI
  FF90 = ifort -lmpi
endif

FF90_FLAGS= -O3

F2PY = f2py
# Note: The F2PY_FF90 variable is ignored for numpy>=2.0 or python>=3.12. FF90 is used instead.
F2PY_FF90 = gnu95

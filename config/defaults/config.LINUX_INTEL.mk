FCOMPILER=mpiifort
FCOMPILER_ALL_FLAGS=$(FCOMPILER) -O3
F2PY_INCLUDES=-I$(MPI_INSTALL_DIR)/include
F2PY_ALL_FLAGS=--fcompiler=intelem --f90exec=$(FCOMPILER) --f77exec=$(FCOMPILER) $(F2PY_INCLUDES) --opt='-O3'
F2PY=python3 -m numpy.f2py

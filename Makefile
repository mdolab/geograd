# makefile for geograd
ARCH_SPECIFIC = config/config.mk

include ${ARCH_SPECIFIC}

# Need to extract MPI include and library flags since Meson for some reason does not link mpi prop
# The following attempts to auto-detect flags for both OpenMPI and Intel MPI
MPI_SHOW := $(shell $(FF90) -show 2>/dev/null)
MPI_INCLUDE_FLAGS := $(shell echo "$(MPI_SHOW)" | grep -oE '(-I[^ ]+)' | tr '\n' ' ')
MPI_LINK_FLAGS := $(shell echo "$(MPI_SHOW)" | grep -oE '(-L[^ ]+|-l[^ ]+|-Wl,[^ ]+)' | tr '\n' ' ')
# The following adds rpath flags if they are not already present
ifeq (,$(findstring -Wl,-rpath,$(MPI_LINK_FLAGS)))
	MPI_RPATH := $(shell echo "$(MPI_LINK_FLAGS)" | grep -oE '\-L[^ ]+' | sed 's/-L/-Wl,-rpath,/' | tr '\n' ' ')
	MPI_LINK_FLAGS += $(MPI_RPATH)
endif

# Include flags for compilation
FF90_INCLUDE_FLAGS = $(MPI_INCLUDE_FLAGS)

# Check what build backend we use for f2py compilation
USE_MESON := $(shell python -c "import sys, numpy; print(int(sys.version_info[:2] >= (3,12) or int(numpy.__version__.split('.')[0]) >= 2))")

F2PY_COMMON := --f90flags="$(FF90_FLAGS) -DINSTRUMENTATION $(FF90_INCLUDE_FLAGS)"
F2PY_COMMON_COMPLEX := --f90flags="$(FF90_FLAGS) -DUSE_COMPLEX $(FF90_INCLUDE_FLAGS)"

ifeq ($(USE_MESON),1)
	export FC=$(FF90)
	export LDFLAGS=$(MPI_LINK_FLAGS)
	F2PY_ALL_FLAGS = $(F2PY_COMMON)
	F2PY_ALL_FLAGS_COMPLEX = $(F2PY_COMMON_COMPLEX)
else
	F2PY_ALL_FLAGS = --fcompiler=$(F2PY_FF90) --f90exec=$(FF90) --f77exec=$(FF90) $(F2PY_COMMON)
	F2PY_ALL_FLAGS_COMPLEX = --fcompiler=$(F2PY_FF90) --f90exec=$(FF90) --f77exec=$(FF90) $(F2PY_COMMON_COMPLEX)
endif


default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_INTEL.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. With the config file specified, rerun "; \
	echo "'make' and the build will start"; \
	else make default_build;\
	fi;

clean:
	find . -name '*.mod' -delete
	find . -name '*.so' -delete
	find . -name '*.msg' -delete
	find . -name '*.msg~' -delete
	find . -name '*.f90~' -delete
	find . -name 'newFile' -delete
	find . -name '*_complex.F90' -delete
	find . -name '*_complex.F90' -delete

test:
	testflo -n 4

default_build: python python_complex

pyf: tapenade/triangles_d.f90 tapenade/triangles_b.f90 src/triangles.F90 src/geograd_parallel.F90 src/geograd.F90
	$(F2PY) $(F2PY_ALL_FLAGS) tapenade/triangles_d.f90 tapenade/triangles_b.f90 src/triangles.F90 src/geograd.F90 src/geograd_parallel.F90 -m geograd -h f2py/geograd.pyf

pyf_complex: complex/triangles_complex.F90 complex/geograd_complex.F90
	$(F2PY) $(F2PY_ALL_FLAGS_COMPLEX) complex/triangles_complex.F90 complex/geograd_complex.F90 -m geograd_complex -h f2py/geograd_complex.pyf

pyf_test: src/mpitest.F90
	$(F2PY) $(F2PY_ALL_FLAGS) src/mpitest.F90 -m mpitest -h f2py/mpitest.pyf

mpitest: src/mpitest.F90 f2py/mpitest.pyf
	$(F2PY) $(F2PY_ALL_FLAGS) -c f2py/mpitest.pyf src/mpitest.F90
	mpirun -np 4 python test_mpi.py

python_complex: complex/triangles_complex.F90 complex/geograd_complex.F90 f2py/geograd_complex.pyf complexify.mod complex/geograd_parallel_complex.F90
	$(F2PY) $(F2PY_ALL_FLAGS_COMPLEX) \
	-c f2py/geograd_complex.pyf \
	complex/complexify.F90 \
	complex/triangles_complex.F90 \
	complex/geograd_complex.F90 complex/geograd_parallel_complex.F90
	mv *.so geograd/libgeograd_complex.so


python: f2py/geograd.pyf src/triangles.F90 src/geograd.F90 triangles_d.mod triangles_b.mod src/geograd_parallel.F90
	$(F2PY) $(F2PY_ALL_FLAGS) \
	-c f2py/geograd.pyf \
	tapenade/triangles_d.f90 tapenade/triangles_b.f90 src/triangles.F90 \
	src/geograd.F90 src/geograd_parallel.F90 \
	tapenade/adBuffer.f tapenade/adStack.c
	mv *.so geograd/libgeograd.so

complex/triangles_complex.F90: src/triangles.F90
	python complex/complexify.py src/triangles.F90
	mv newFile complex/triangles_complex.F90

complex/geograd_complex.F90: src/geograd.F90
	python complex/complexify.py src/geograd.F90
	mv newFile complex/geograd_complex.F90

complex/geograd_parallel_complex.F90: src/geograd_parallel.F90
	python complex/complexify.py src/geograd_parallel.F90
	mv newFile complex/geograd_parallel_complex.F90

complexify.mod: complex/complexify.F90
	$(FF90) $(FF90_FLAGS) -c complex/complexify.F90
	rm complexify.o

triangles_d.mod: tapenade/triangles_d.f90
	$(FF90) $(FF90_FLAGS) -c tapenade/triangles_d.f90
	rm triangles_d.o

triangles_b.mod: tapenade/triangles_b.f90
	$(FF90) $(FF90_FLAGS) -c tapenade/triangles_b.f90
	rm triangles_b.o

tapenade: src/triangles.F90
	tapenade src/triangles.F90 -d -root point_tri -root line_line -root intersect
	tapenade src/triangles.F90 -b -root point_tri -root line_line -root intersect
	mv triangles_d.f90 tapenade/triangles_d.f90
	mv triangles_b.f90 tapenade/triangles_b.f90

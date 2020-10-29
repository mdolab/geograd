# makefile for geograd
ARCH_SPECIFIC = config/config.mk

include ${ARCH_SPECIFIC}

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
	python3 tests/test_primitives.py
	python3 tests/test_integration.py
	mpirun -np 4 python3 tests/test_integration_parallel.py

default_build: python3 python3_complex test

pyf: tapenade/triangles_db.f90 src/triangles.F90 src/geograd_parallel.F90 src/geograd.F90
	$(F2PY) tapenade/triangles_db.f90 src/triangles.F90 src/geograd.F90 src/geograd_parallel.F90 -m geograd -h f2py/geograd.pyf

pyf_complex: complex/triangles_complex.F90 complex/geograd_complex.F90
	$(F2PY) complex/triangles_complex.F90 complex/geograd_complex.F90 -m geograd_complex -h f2py/geograd_complex.pyf

pyf_test: src/mpitest.F90
	$(F2PY) src/mpitest.F90 -m mpitest -h f2py/mpitest.pyf

mpitest: src/mpitest.F90 f2py/mpitest.pyf
	$(F2PY) $(F2PY_ALL_FLAGS) -c f2py/mpitest.pyf src/mpitest.F90
	mpirun -np 4 python3 test_mpi.py

python3_complex: complex/triangles_complex.F90 complex/geograd_complex.F90 f2py/geograd_complex.pyf complexify.mod complex/geograd_parallel_complex.F90
	$(F2PY) $(F2PY_ALL_FLAGS) -DUSE_COMPLEX -c f2py/geograd_complex.pyf complex/triangles_complex.F90 complex/geograd_complex.F90 complex/complexify.F90 complex/geograd_parallel_complex.F90
	mv *.so geograd/libgeograd_complex.so

python3: tapenade/triangles_db.f90 f2py/geograd.pyf src/triangles.F90 src/geograd.F90 triangles_db.mod src/geograd_parallel.F90
	$(F2PY) $(F2PY_ALL_FLAGS) -DINSTRUMENTATION -c f2py/geograd.pyf tapenade/triangles_db.f90 src/triangles.F90 src/geograd.F90 src/geograd_parallel.F90 tapenade/adBuffer.f tapenade/adStack.c
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
	$(FCOMPILER_ALL_FLAGS) -c complex/complexify.F90
	rm complexify.o

triangles_db.mod: tapenade/triangles_db.f90
	$(FCOMPILER_ALL_FLAGS) -c tapenade/triangles_db.f90
	rm triangles_db.o

tapenade/triangles_db.f90: src/triangles.F90
	tapenade src/triangles.F90 -d -b -root point_tri -root line_line -root intersect
	mv triangles_db.f90 tapenade/triangles_db.f90

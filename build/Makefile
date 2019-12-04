FCOMPILER=mpifort
F2PY_INCLUDES=-I$(MPI_INSTALL_DIR)/include
F2PY_ALL_FLAGS=--f90exec=$(FCOMPILER) $(F2PY_INCLUDES)
F2PY=python3 -m numpy.f2py 

default: python3 python3_complex test

pyf: tapenade/triangles_db.f90 src/triangles.F90
	$(F2PY) tapenade/triangles_db.f90 src/triangles.F90 src/geograd.F90 -m geograd -h f2py/geograd.pyf

pyf_complex: complex/triangles_complex.F90 complex/geograd_complex.F90
	$(F2PY) complex/triangles_complex.F90 complex/geograd_complex.F90 -m geograd_complex -h f2py/geograd_complex.pyf

pyf_test: src/mpitest.F90
	$(F2PY) src/mpitest.F90 -m mpitest -h f2py/mpitest.pyf

mpitest: src/mpitest.F90 f2py/mpitest.pyf
	$(F2PY) $(F2PY_ALL_FLAGS) -c f2py/mpitest.pyf src/mpitest.F90
	mpirun -np 4 python3 test_mpi.py

python3_complex: complex/triangles_complex.F90 complex/geograd_complex.F90 f2py/geograd_complex.pyf complexify.mod
	$(F2PY) $(F2PY_ALL_FLAGS) -DUSE_COMPLEX -c f2py/geograd_complex.pyf complex/triangles_complex.F90 complex/geograd_complex.F90 complex/complexify.F90

python3: tapenade/triangles_db.f90 f2py/geograd.pyf src/triangles.F90 src/geograd.F90 triangles_db.mod
	$(F2PY) $(F2PY_ALL_FLAGS) -c f2py/geograd.pyf tapenade/triangles_db.f90 src/triangles.F90 src/geograd.F90 tapenade/adBuffer.f tapenade/adStack.c

complex/triangles_complex.F90: src/triangles.F90
	python complex/complexify.py src/triangles.F90
	mv newFile complex/triangles_complex.F90

complex/geograd_complex.F90: src/geograd.F90
	python complex/complexify.py src/geograd.F90
	mv newFile complex/geograd_complex.F90

complexify.mod: complex/complexify.F90
	gfortran -c complex/complexify.F90
	rm complexify.o

triangles_db.mod: tapenade/triangles_db.f90
	gfortran -c tapenade/triangles_db.f90
	rm triangles_db.o

tapenade/triangles_db.f90: src/triangles.F90
	tapenade src/triangles.F90 -d -b -root point_tri -root line_line -root intersect
	mv triangles_db.f90 tapenade/triangles_db.f90

clean:
	find . -name '*.mod' -delete
	find . -name '*.so' -delete 
	find . -name '*.msg' -delete
	find . -name '*.msg~' -delete
	find . -name '*.f90~' -delete
	find . -name 'newFile' -delete
	find . -name '*_db.f90' -delete 
	find . -name '*_complex.F90' -delete
	find . -name '*_complex.F90' -delete

test:
	python3 test_primitives.py
	python3 test_integration.py
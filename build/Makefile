default: python3 python3_complex test

pyf: tapenade/triangles_db.f90 src/triangles.F90
	python3 -m numpy.f2py tapenade/triangles_db.f90 src/triangles.F90 src/geograd.F90 -m geograd -h f2py/geograd.pyf

pyf_complex: complex/triangles_complex.F90 complex/geograd_complex.F90
	python3 -m numpy.f2py complex/triangles_complex.F90 complex/geograd_complex.F90 -m geograd_complex -h f2py/geograd_complex.pyf

python3_complex: complex/triangles_complex.F90 complex/geograd_complex.F90 f2py/geograd_complex.pyf complex/complexify.mod
	python3 -m numpy.f2py -c f2py/geograd_complex.pyf complex/triangles_complex.F90 complex/geograd_complex.F90 complex/complexify.F90

python3: tapenade/triangles_db.f90 f2py/geograd.pyf src/triangles.F90 src/geograd.F90
	python3 -m numpy.f2py -c f2py/geograd.pyf tapenade/triangles_db.f90 src/triangles.F90 src/geograd.F90 tapenade/adBuffer.f tapenade/adStack.c

complex/triangles_complex.F90: src/triangles.F90
	python complex/complexify.py src/triangles.F90
	mv newFile complex/triangles_complex.F90

complex/geograd_complex.F90: src/geograd.F90
	python complex/complexify.py src/geograd.F90
	mv newFile complex/geograd_complex.F90

complex/complexify.mod: complex/complexify.F90
	gfortran -c complex/complexify.F90
	rm complexify.o

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
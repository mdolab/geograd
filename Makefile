


default: triangles_db.f90 triangles.pyf triangles.F90
	python3 -m numpy.f2py -c triangles.pyf triangles_db.f90 triangles.F90 adBuffer.f adStack.c

complexify: triangles.F90
	python complexify.py triangles.F90
	mv newFile triangles_complex.F90

pyf_complex: triangles_complex.F90
	python3 -m numpy.f2py triangles_complex.F90 -m triangles_complex -h triangles_complex.pyf

python3_complex: triangles_complex.F90 triangles_complex.pyf 
	python3 -m numpy.f2py -c triangles_complex.pyf triangles_complex.F90 complexify.F90

python3: triangles_db.f90 triangles.pyf triangles.F90
	python3 -m numpy.f2py -c triangles.pyf triangles_db.f90 triangles.F90 adBuffer.f adStack.c

pyf: triangles_db.f90
	python3 -m numpy.f2py triangles_db.f90 triangles.F90 -m triangles -h triangles.pyf

tapenade: triangles.F90
	tapenade triangles.F90 -d -b -root point_tri -root line_line

clean:
	rm *.mod -f
	rm *.so -f 
	rm *.msg -f
	rm *.msg~ -f
	rm *.f90~ -f
	rm newFile -f
	rm *_db.f90 -f 
	rm *_complex.F90 -f

test:
	python3 test.py
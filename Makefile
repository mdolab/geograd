


default: triangles.F90
	python -m numpy.f2py -c triangles.F90 -m triangles

python3: triangles_db.f90 triangles.pyf
	python3 -m numpy.f2py -c triangles.pyf triangles_db.f90 triangles.F90

pyf: triangles_db.f90
	python3 -m numpy.f2py triangles_db.f90 triangles.F90 -m triangles -h triangles.pyf

tapenade: triangles.F90
	tapenade triangles.F90 -d -b -root point_tri -root line_line

clean:
	rm *.mod -f
	rm *.so -f 
	rm *.msg~ -f
	rm *.f90~ -f

test:
	python3 test.py



default: triangles.F90
	python -m numpy.f2py -c triangles.F90 -m triangles

python3: triangles.F90
	python3 -m numpy.f2py -c triangles.F90 -m triangles

clean:
	rm *.mod -f
	rm *.so -f 

test:
	python test.py
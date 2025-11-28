# Include user supplied configuration file
include ../../config/config.mk

# Group Fortran compiler flags
FF90_ALL_FLAGS = $(FF90_FLAGS) -DUSE_COMPLEX
CC_ALL_FLAGS = $(C_FLAGS)
LINKER_ALL_FLAGS = $(FF90_FLAGS) $(LINKER_FLAGS)


# Include full list of files
include fileList

# Include full list of directories
include directoryList

# Include the list of rules
include ../../src/build/rules

# Prepend complexify.F90 to the fortranFiles list to ensure it is built first
fortranFiles := complexify.F90 $(fortranFiles)

# Set the make VPATH variable to the "dirs" variable from
# directorylist. We must first append the '../'
dirs:=$(addprefix ../,$(dirs))
VPATH:=$(dirs)

# Add ../ prefix to all source files, remove directories
fortranFiles:=$(addprefix ../,$(fortranFiles))
cFiles:=$(addprefix ../,$(cFiles))

fortranFilesNoDir=$(notdir $(fortranFiles))
cFilesNoDir=$(notdir $(cFiles))

# Generate separate list of .F90 and .f90 files using the filter command
f90Files=$(filter %.f90,$(fortranFilesNoDir))
F90Files=$(filter %.F90,$(fortranFilesNoDir))

# Finally convert all source files to .o
OFILES=$(f90Files:%.f90=%.o) $(F90Files:%.F90=%.o) $(cFilesNoDir:%.c=%.o)

.PHONY: sources lib

PYTHON_OBJECTS = fortranobject.o \
		 libgeograd_complexmodule.o\
		 libgeograd_complex-f2pywrappers2.o\


default: lib ../f2py/geograd_complex.pyf

# Generate Python inlude directory
	$(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print(numpy.get_include())'))
	 @echo "#------------------------------------------------------#"
	 @echo Numpy Include Directory: $(NUMPY_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) ../../src/f2py/get_f2py.py))
	 @echo "#------------------------------------------------------#"
	 @echo f2py root directory: $(F2PY_ROOT)
	 @echo "#------------------------------------------------------#"

# Run f2py
	$(F2PY) ../f2py/geograd_complex.pyf

# Compile f2py-generated C wrapper file
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libgeograd_complexmodule.c

# Compile the Fortran wrapper in the f2py source
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) -c \
	$(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compile f2py-generated Fortran wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I./ -c libgeograd_complex-f2pywrappers2.f90

# Create shared library
	$(FF90) -shared $(PYTHON_OBJECTS) libgeograd_complex.a $(LINKER_ALL_FLAGS) -o libgeograd_complex.so

# Test if the shared library can be imported in Python
	$(PYTHON) ../../src/f2py/importTest.py libgeograd_complex.so

# Move the shared library to the main Python directory
	mv libgeograd_complex.so ../../geograd

sources: $(OFILES)

lib: sources
	$(AR) $(AR_FLAGS) libgeograd_complex.a $(OFILES)

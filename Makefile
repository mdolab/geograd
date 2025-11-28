# makefile for geograd

default:
# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. With the config file specified, rerun "; \
	echo "'make' and the build will start"; \
	else make geograd_build;\
	fi;

clean:
	rm -fr src/build/*.mod
	rm -fr src/build/*.o
	rm -fr src/build/*.a
	rm -fr src/build/*.so
	rm -f *~ config.mk;

test:
	testflo -n 4

geograd_build:
	ln -sf config/config.mk config.mk;
	(cd src/build/ && make)
	(cd src_cs/build/ && make)

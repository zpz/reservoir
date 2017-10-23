# Do not add machine specific lib paths.
#
# Instead of including required paths in the makefile,
# add them to the environ variable
#
#   DYLD_FALLBACK_LIBRARY_PATH
#
# on a Mac and
#
#   LD_LIBRARY_PATH
#
# on Linux.


CC = g++-4.7
#CCFLAGS := -std=c++11 -Wfatal-errors -ggdb3 -fPIC -Wall
CCFLAGS = -std=c++11 -Wfatal-errors -ggdb3 -fPIC -Wall -O2
LLFLAGS = -m64 -shared
RES_LIBS = -L$(HOME)/usr/lib -lhdf5util -lhdf5_hl -lhdf5
RES_INCLUDES = -I$(HOME)/usr/include
H5_LIBS = -lhdf5_hl -lhdf5
H5_INCLUDES =
INSTALLDIR = $(HOME)/usr


all: libhdf5util.so libreservoir.so


libhdf5util.so: hdf5util.o
	$(CC) $(LLFLAGS) $^ $(H5_LIBS) -o $@
	install $@ $(INSTALLDIR)/lib/
	cp -f hdf5util.h $(INSTALLDIR)/include/

hdf5util.o: hdf5util.cpp hdf5util.h
	$(CC) $(CCFLAGS) $(H5_INCLUDES) -c $< -o $@

libreservoir.so: reservoir.o
	$(CC) $(LLFLAGS) -o $@ $^ $(RES_LIBS)
	install $@ $(INSTALLDIR)/lib/
	cp -f reservoir.h $(INSTALLDIR)/include/

reservoir.o: reservoir.cpp reservoir.h
	$(CC) $(CCFLAGS) $(RES_INCLUDES) -c $< -o $@

clean:
	rm -f hdf5util.o reservoir.o
	rm -f libhdf5util.so libreservoir.so
	rm -f $(INSTALLDIR)/lib/libhdf5util.so
	rm -f $(INSTALLDIR)/include/hdf5util.h
	rm -f $(INSTALLDIR)/lib/libreservoir.so
	rm -f $(INSTALLDIR)/include/reservoir.h


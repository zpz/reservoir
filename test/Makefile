CC = g++-4.7
#CCFLAGS := -std=c++11 -Wfatal-errors -ggdb3 -Wall
CCFLAGS = -std=c++11 -Wfatal-errors -ggdb3 -Wall -O2
LLFLAGS =
RES_INCLUDES = -I../ -I./ -I$(HOME)/usr/include
RES_LIBS = -L$(HOME)/usr/lib -lrandom -lhdf5util
H5_INCLUDES = -I../
H5_LIBS = -L../ -lhdf5util -lhdf5_hl -lhdf5

all: test_reservoir test_h5

test_reservoir: test_reservoir.o
	$(CC) $(LLFLAGS) $^ $(RES_LIBS) -o $@

test_reservoir.o: %.cpp ../reservoir.h ../libreservoir.so
	$(CC) -c $(CCFLAGS) $(RES_INCLUDES) $< -o $@

test_h5: test_h5.o
	$(CC) $(LLFLAGS) $^ $(H5_LIBS) -o $@

test_h5.o: %.cpp ../hdf5util.h ../libhdf5util.so
	$(CC) -c $(CCFLAGS) $(H5_INCLUDES) $< -o $@


clean:
	rm -f *.o
	rm -f test_reservoir test_h5
	rm -f *h5


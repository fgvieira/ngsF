CC=gcc
CXX=g++

#CFLAGS = -g -Wall
CFLAGS = -O3 -Wall
DFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE #-D_USE_BGZF
LIB = -lgsl -lgslcblas -lz -lpthread

all: bgzip ngsF


bgzip:
	$(MAKE) -C bgzf bgzip



parse_args.o: parse_args.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

read_data.o: read_data.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c read_data.cpp

EM.o: EM.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c EM.cpp

shared.o: shared.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c shared.cpp

ngsF: ngsF.cpp parse_args.o read_data.o EM.o shared.o
	$(CXX) $(CFLAGS) $(DFLAGS) ngsF.cpp parse_args.o read_data.o EM.o shared.o bgzf/bgzf.o bgzf/knetfile.o $(LIB) -o ngsF

test:
	@cd examples/; bash test.sh 2> test.log; cd ../

clean:
	@rm -f bgzf/*.o bgzf/bgzip *.o ngsF examples/testF.*


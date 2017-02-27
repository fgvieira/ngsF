CC?=gcc
CXX?=g++

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
all: ngsF
else
HTS_INCDIR=$(realpath ./htslib)
HTS_LIBDIR=$(realpath ./htslib/bgzf.o) $(realpath ./htslib/knetfile.o)
all: bgzip ngsF
endif

#CFLAGS = -g -Wall
CFLAGS = -O3 -Wall
DFLAGS = -I$(HTS_INCDIR) -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE #-D_USE_BGZF
LIB = $(shell pkg-config --cflags --libs gsl) -lz -lpthread $(HTS_LIBDIR)



bgzip:
	$(MAKE) -C htslib bgzip



parse_args.o: parse_args.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

read_data.o: read_data.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c read_data.cpp

EM.o: EM.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c EM.cpp

shared.o: shared.cpp shared.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c shared.cpp

ngsF: ngsF.cpp parse_args.o read_data.o EM.o shared.o
	$(CXX) $(CFLAGS) $(DFLAGS) ngsF.cpp parse_args.o read_data.o EM.o shared.o $(LIB) -o ngsF

test:
	@cd examples/; bash test.sh 2> test.log; cd ../

clean:
	@rm -f htslib/*.o htslib/bgzip *.o ngsF examples/testF.*


INCLUDES=
LIBS=
CXXFLAGS=-Wall -pedantic -std=c++11 -O3 $(INCLUDES)
LDFLAGS=$(LIBS) -lsdsl -ldivsufsort -ldivsufsort64 -ljellyfish-2.0 \
		-lpthread -lz

CXX=g++

all: lean df

sbf: main.o Kmers.o SBF.o SBF1.o SBF2.o BF.o util.o Count.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o sbf *~
# DO NOT DELETE

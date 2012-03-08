
EXECUTABLES = illuminus
TARGETS = $(EXECUTABLES)

RNG := $(patsubst %.c,%.o,$(wildcard other_libraries/rng/*.c))
NM := $(patsubst %.cpp,%.o,$(wildcard other_libraries/newmat11/*.cpp))

AR = ar
CXX = g++
CXXFLAGS = -O3 -Wno-deprecated -I./
LIBPATH = -L./other_libraries/rng/ -L./other_libraries/newmat11/ -L./
LDFLAGS = $(LIBPATH) -lm -lnewmat -lrng -lplinkbin

.PHONY: test

all: $(EXECUTABLES)

illuminus: illuminus.o librng.a libnewmat.a
	$(CXX) -static $< $(LDFLAGS) -o $@

illuminus.o: illuminus.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

librng.a : $(RNG)
	$(AR) rc $@ $^

libnewmat.a : $(NM)
	$(AR) rc $@ $^

clean : 
	rm -f *.o ./other_libraries/rng/*.o ./other_libraries/newmat11/*.o *.a
	rm -f $(EXECUTABLES)

LIBS = -L./other_libraries/rng/ -L./other_libraries/newmat11/ -L/software/varinf/gftools/lib -lplinkbin
CFLAGS = -O3
CXXFLAGS = -O3 -Wno-deprecated

CXX = /usr/bin/g++

RNG := $(patsubst %.c,%.o,$(wildcard other_libraries/rng/*.c))
NM := $(patsubst %.cpp,%.o,$(wildcard other_libraries/newmat11/*.cpp))

illuminus: illuminus.o librng libnewmat 
	$(CXX) illuminus.o  -o illuminus -lm -lstdc++ -lnewmat -lrng $(LIBS)

illuminus.o: illuminus.cc illuminus.h
	/usr/bin/g++ -O3 -Wno-deprecated -I/software/varinf/gftools/include -c -o illuminus.o illuminus.cc

librng : $(RNG) 
	/usr/bin/ar rc other_libraries/rng/librng.a $(RNG)

libnewmat : $(NM) 
	/usr/bin/ar rc other_libraries/newmat11/libnewmat.a $(NM)

clean : 
	rm -f *.o ./other_libraries/rng/*.o ./other_libraries/rng/*.a ./other_libraries/newmat11/*.o ./other_libraries/newmat11/*.a; rm illuminus
				

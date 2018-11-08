#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.
#O2?

MODE=-Ofast
STD=-std=c++11
DESTDIR=/usr/local/bin
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.0\"

MKLROOT=/software/intel-tools-2015/mkl

INCLUDE=-I ./ -I ../common/include/ -I ~/src/lz4/lib  -I ~/src/eigen-eigen-5a0156e40feb #-I${MKLROOT}/include 
								        #|#
#LDLIBS=-lgsl -lblas -L ../common -lmapgd -L ~/src/lz4/lib -llz4 -fopenmp -Bstatic -lmkl_intel_lp64 -lmkl_core -liomp5 -lmkl_intel_thread -lpthread -lm -ldl -m64 -L${MKLROOT}/lib/intel64/ -L/software/intel-tools-2015/composer_xe_2015.0.090/compiler/lib/intel64/ #-static
LDLIBS=-lgsl -lblas -L ../common -lmapgd -L ~/src/lz4/lib -llz4 -fopenmp -lpthread -lm -ldl #-static

DEBUG=-g -D DEBUG

CXXFLAGS=$(VERSION) $(STD) $(WARNINGS) -DEIGEN_USE_MKL_ALL -DNDEBUG $(INCLUDE) -fopenmp

#CXXFLAGS=$(VERSION) $(STD) $(WARNINGS) -DNDEBUG -frounding-math -pedantic $(INCLUDE) -fopenmp #$(DEBUG) #-mpopcnt
CXXFLAGS=$(VERSION) $(STD) $(WARNINGS) -DNDEBUG -pedantic $(INCLUDE) -fopenmp #$(DEBUG) #-mpopcnt

#CXXFLAGS=$(VERSION) $(STD) $(DEBUG) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) -fopenmp

OBJDIR=obj

TARGET=doge

VPATH = ./src

ifeq ($(HAVE_MAN), yes)
	CXXFLAGS += -D USE_MAN -D PATH=$(DOCPATH)
endif

O2=pseudo.o \
	sketch.o \
	ll.o

#Ofast=mle_5.o
#Ofast=mle_5_scan.o
Ofast=mle_6.o
#Ofast=lagrange_mle.o

O2wDIR=$(addprefix $(OBJDIR)/, $(O2) )
OfastwDIR=$(addprefix $(OBJDIR)/, $(Ofast) )

OBJ=$(O2wDIR) $(OfastwDIR)

#>! Just a list of flags to try at some later date.
# -Ofast -ffloat-store -fexcess-precision=style -ffast-math -fno-rounding-math -fno-signaling-nans
# -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math
# -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math
# -frounding-math -fsingle-precision-constant -fcx-limited-range -fcx-fortran-rules

#PGO

#GCC has Profile-Guided Optimisations features. There isn't a lot of precise GCC documentation about it, but nevertheless getting it to run is quite straightforward.

#first compile your program with -fprofile-generate.
#let the program run (the execution time will be significantly slower as the code is also generating profile information into .gcda files).
#recompile the program with -fprofile-use. If your application is multi-threaded also add the -fprofile-correction flag

#>! TODO ORGANIZE THESE!

all : $(TARGET)

./obj/pseudo.o : pseudo.cc
	$(CXX) -O2 $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

./obj/ll.o : ll.cc
	$(CXX) -O2 $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

./obj/sketch.o : sketch.cc
	$(CXX) -O2 $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

./obj/mle_6.o : mle_6.cc
	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

./obj/lagrange_mle.o : lagrange_mle.cc
	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

./obj/mle_5.o : mle_5.cc
	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

./obj/mle_5_scan.o : mle_5_scan.cc
	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)
	
#$(O2wDIR) : $(O2:.o=.cc)
#	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

#$(OfastwDIR) : $(Ofast:.o=.cc)
#	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

#./obj/%.o : %.cc
#	$(CXX) -Ofast $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

$(TARGET) : $(OBJ)
	$(CXX) -Ofast $(CXXFLAGS) $(OBJ) -o $(TARGET) $(LDLIBS) 

clean:
	rm -f $(OBJDIR)/*.o

.PHONY: clean

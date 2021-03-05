# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O
#     parallelization by OpenMP but no MPI dependencies
# Change these lines to fit your system.
CPP = g++
MPICPP = mpic++

fortranint=4
BLASLIB = /usr/lib/x86_64-linux-gnu/libblas.so.3
LAPACKLIB = /usr/lib/x86_64-linux-gnu/liblapack.so.3



# The rest should not have to be changed.

UTILLIB = ../util/lib/libutil_$(mode).a
UTILINC = -I ../util/include

INCLUDE = -I src -I include -I .. -I . $(UTILINC)
LIBS = $(LAPACKLIB) $(UTILLIB) $(BLASLIB) -lgfortran
LIBDIRS =
ifneq (,$(findstring O,$(opt)))
	OPTFLAGS = $(opt) #-DNDEBUG
else
	OPTFLAGS = $(opt)
endif
ifeq ($(mode),s)
	OMP = -fopenmp -DUTIL_OMP
else
	OMP =
endif
CPPFLAGS = $(OPTFLAGS) $(OMP) -DFORTRAN_INT_$(fortranint) -DAMB_SPECIAL #-DANALYZE_CODE #-Wall #-Wextra
LDFLAGS = $(OMP)

.SUFFIXES:
.SUFFIXES: .cpp .o

CCPPSRCS = src/Hd.cpp src/Compress.cpp src/Hmat.cpp src/HmatIo.cpp

COBJECTS = $(patsubst %.cpp,%.o,$(CCPPSRCS))

%.o : %.cpp
ifeq ($(mode),p)
	$(MPICPP) -DUTIL_MPI $(CPPFLAGS) $(INCLUDE) -c $< -o $@
else
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
endif

all: libhmmvp

libhmmvp: $(COBJECTS)
	ar rucs lib/libhmmvp_$(mode).a $(COBJECTS)

build: libhmmvp
ifeq ($(mode),s)
	$(CPP) src/hmmvpbuild.cpp $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o bin/hmmvpbuild
else
	$(MPICPP) src/hmmvpbuild.cpp $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) -DUTIL_MPI $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o bin/hmmvpbuild
endif

test: libhmmvp
	$(CPP) check_mvp_time_scaling.cpp $(CPPFLAGS) $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o check_mvp_time_scaling
	$(CPP) test_zero.cpp $(CPPFLAGS) $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o test_zero

clang :
	clang -fsyntax-only $(INCLUDE) src/*.cpp

clean:
	rm -f src/*.o external/*.o

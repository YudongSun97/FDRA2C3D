# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O
#     serial with no MPI dependencies
# Change these lines to fit your system.
BLASLIB = /cm/shared/apps/blas/open64/1/lib64/libblas.a
LAPACKLIB = /cm/shared/apps/lapack/open64/64/3.5.0/liblapack.so
#BLASLIB = /usr/lib/x86_64-linux-gnu/libblas.so.3
#LAPACKLIB = /usr/lib/x86_64-linux-gnu/liblapack.so.3
#BLASLIB = ~/dl/lapack-3.2.1/libblas.a
#LAPACKLIB = ~/dl/lapack-3.2.1/liblapack.a
CPP = g++
MPICPP = mpic++
FORTRAN = gfortran
# The rest should not have to be changed.

INCLUDE = -I src -I include -I ..
LIBS = $(LAPACKLIB) $(BLASLIB) -lgfortran
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
CPPFLAGS = $(OPTFLAGS) $(OMP) #-Wall #-Wextra #-DANALYZE_CODE
FFLAGS = $(OPTFLAGS) $(OMP) -Wall
LDFLAGS = $(OMP)

.SUFFIXES:
.SUFFIXES: .cpp .f .o

CCPPSRCS =  src/Mpi.cpp src/KeyValueFile.cpp src/CodeAnalysis.cpp \
src/ValueSetter.cpp src/PolyInterp.cpp
CFSRCS = 

COBJECTS = $(patsubst %.cpp,%.o,$(CCPPSRCS)) $(patsubst %.f,%.o,$(CFSRCS))

%.o : %.cpp
ifeq ($(strip $(mode)),p)
	$(MPICPP) -DUTIL_MPI $(CPPFLAGS) $(INCLUDE) -c $< -o $@
else
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
endif

%.o : %.f
	$(FORTRAN) $(FFLAGS) -c $< -o $@

all : libutil

libutil : $(COBJECTS)
	ar rucs lib/libutil_$(mode).a $(COBJECTS)

clang :
	clang -fsyntax-only $(INCLUDE) src/*.cpp dc3dm/*.cpp

clean:
	rm -f src/*.o dc3dm/*.o external/*.o

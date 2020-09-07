# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O
#     serial with no MPI dependencies
# Change these lines to fit your system.
HMMVPINC = -I /home/camcat/code/ambrad
HMMVPLIB = /home/camcat/code/hmmvp/lib/libhmmvp_$(mode).a
UTILINC = -I /home/camcat/code/ambrad
UTILLIB = /home/camcat/code/lib/libutil_$(mode).a
BLASLIB = /home/camcat/code/lib/lapack-3.2.1/libblas.a

CPP = g++
MPICPP = mpic++
# The rest should not have to be changed.

CPPSRCS = src/main.cpp src/OdeFn.cpp src/Ode.cpp src/StreamDataFile.cpp	\
src/Fdra.cpp src/HmatrixStressFn.cpp

INCLUDE = -I src -I include $(HMMVPINC) $(UTILINC)
LIBS = $(HMMVPLIB) $(UTILLIB) $(BLASLIB) -lgfortran -lgomp
LIBDIRS =
OPTFLAGS = $(opt)
CPPFLAGS = $(OPTFLAGS) -Wall #-DANALYZE_CODE
LDFLAGS = $(OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJECTS = $(patsubst %.cpp,%.o,$(CPPSRCS))

%.o: %.cpp
ifeq ($(mode),s)
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
else
	$(MPICPP) -DUTIL_MPI $(CPPFLAGS) $(INCLUDE) -c $< -o $@
endif

all: $(OBJECTS)
ifeq ($(mode),s)
	$(CPP) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) $(OBJECTS) $(LIBS) -o bin/fdra
else
	$(MPICPP) $(LDFLAGS) $(LIBFLAGS) -DUTIL_MPI $(LIBDIRS) $(OBJECTS) $(LIBS) -o bin/fdra
endif

add_bc:
	$(CPP) test/add_bc_to_fdra_kvf.cpp $(UTILLIB) $(INCLUDE) -o bin/add_bc_to_fdra_kvf

fdraclang:
	clang -fsyntax-only $(INCLUDE) $(CPPSRCS)

clean:
	rm -f src/*.o hmmvp/*.o bin/fdra

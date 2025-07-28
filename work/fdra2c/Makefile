# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O [debug=1]
#     serial with no MPI dependencies
# Change these lines to fit your system.


# HMMVP has one serial and two parallel modes: "omp" and "mpi"
ifeq ($(mode),p)
  mode_hmmvp = omp
else
  mode_hmmvp = s
endif

HMMVPINC = -I external/hmmvp
HMMVPLIB = external/hmmvp/lib/libhmmvp_$(mode_hmmvp).a
# BLASLIB = /home/camcat/code/lib/lapack-3.2.1/libblas.a
# Yudong Mar 7 2025
BLASLIB = -lopenblas

ifeq ($(debug),1)
  db = "_db"
  DBFLAGS= -g -O0
else
  db = 
  DBFLAGS = 
endif

CPP = g++
MPICPP = mpic++
# The rest should not have to be changed.

CPPSRCS = src/main.cpp src/OdeFn.cpp src/Ode.cpp src/StreamDataFile.cpp	\
src/Fdra.cpp src/HmatrixStressFn.cpp \
external/hmmvp/src/Mpi.cpp external/hmmvp/src/ValueSetter.cpp

CPPSRCS_BC = test/add_bc_to_fdra_kvf.cpp external/hmmvp/src/KeyValueFile.cpp

INCLUDE = -I src -I include $(HMMVPINC)
LIBS = $(HMMVPLIB) $(BLASLIB) -lgfortran -lgomp
LIBDIRS =
OPTFLAGS = $(DBFLAGS) $(opt)
CPPFLAGS = $(OPTFLAGS) #-DANALYZE_CODE
LDFLAGS = $(OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJECTS = $(patsubst %.cpp,%.o,$(CPPSRCS))

%.o: %.cpp

ifeq ($(mode),p)
	$(MPICPP) -DUTIL_MPI $(CPPFLAGS) $(INCLUDE) -c $< -o $@
else
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
endif

all: $(OBJECTS)
ifeq ($(mode),s)
	$(CPP) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) $(OBJECTS) $(LIBS) -o bin/fdra$(db)
else
	$(MPICPP) $(LDFLAGS) $(LIBFLAGS) -DUTIL_MPI $(LIBDIRS) $(OBJECTS) $(LIBS) -o bin/fdra$(db)
endif

add_bc:
	$(CPP) $(CPPFLAGS) $(DBFLAGS) $(CPPSRCS_BC) $(INCLUDE) -o bin/add_bc_to_fdra_kvf$(db)

fdraclang:
	clang -fsyntax-only $(INCLUDE) $(CPPSRCS)

clean:
	rm -f src/*.o hmmvp/*.o bin/fdra$(db)






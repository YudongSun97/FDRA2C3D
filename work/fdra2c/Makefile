# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O [debug=1]
#     serial with no MPI dependencies
# Change these lines to fit your system.

HMMVPINC = -I external
HMMVPLIB = external/hmmvp/lib/libhmmvp_$(mode).a
UTILINC = #-I external/util/include
UTILLIB = external/util/lib/libutil_$(mode).a
BLASLIB = /home/camcat/code/lib/lapack-3.2.1/libblas.a
#BLASLIB = /cm/shared/apps/blas/open64/1/lib64/libblas.a

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
src/Fdra.cpp src/HmatrixStressFn.cpp src/UniLoadStressFn.cpp#\
/home/camcat/code/ambrad/util/src/ValueSetter.cpp \
/home/camcat/code/ambrad/util/src/KeyValueFile.cpp \
/home/camcat/code/ambrad/util/src/Mpi.cpp \
/home/camcat/code/ambrad/util/src/CodeAnalysis.cpp \

INCLUDE = -I src -I include $(HMMVPINC) $(UTILINC)
LIBS = $(HMMVPLIB) $(UTILLIB) $(BLASLIB) -lgfortran -lgomp
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
	$(CPP) $(CPPFLAGS) $(DBFLAGS) test/add_bc_to_fdra_kvf.cpp $(UTILLIB) $(INCLUDE) -o bin/add_bc_to_fdra_kvf$(db)

fdraclang:
	clang -fsyntax-only $(INCLUDE) $(CPPSRCS)

clean:
	rm -f src/*.o hmmvp/*.o bin/fdra$(db)






# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O
#     serial with no MPI dependencies
# Change these lines to fit your system.
CPP = g++
MPICPP = mpic++
# The rest should not have to be changed.

INCLUDE = -I src -I include -I ..
OPTFLAGS = $(opt)
CPPFLAGS = $(OPTFLAGS) -Wall -DANALYZE_CODE
LDFLAGS = $(OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .cpp .o

CCPPSRCS = src/Mpi.cpp src/KeyValueFile.cpp src/CodeAnalysis.cpp	\
src/ValueSetter.cpp src/RectMeshUnstruct.cpp

COBJECTS = $(patsubst %.cpp,%.o,$(CCPPSRCS))

%.o : %.cpp
ifeq ($(mode),p)
	$(MPICPP) -DUTIL_MPI $(CPPFLAGS) $(INCLUDE) -c $< -o $@
else
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
endif

all : libutil meshrect

libutil : $(COBJECTS)
	ar rucs lib/libutil.a $(COBJECTS)

meshrect : libutil
	$(CPP) $(INCLUDE) $(CPPFLAGS) $(LDFLAGS) $(LIBFLAGS) src/MakeRectUnstructMesh.cpp lib/libutil.a -o ./bin/MakeRectUnstructMesh

clang :
	clang -fsyntax-only $(INCLUDE) src/*.cpp

clean:
	rm -f src/*.o lib/libutil.a bin/*

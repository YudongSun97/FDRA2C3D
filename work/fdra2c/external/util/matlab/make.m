% Configure options
use_omp = 1;

p = '.';
flags = '-O';
if(use_omp)
  flags = [flags [' CXXFLAGS="\$CXXFLAGS -fopenmp" ',...
                  'LDFLAGS="\$LDFLAGS -fopenmp"']];
end
flags = [flags ' CXXLIBS="\$CXXLIBS -Wall -lmwblas -lmwlapack"'];
mc = sprintf('mex -I../.. -outdir %s -largeArrayDims %s -DUSING_UNIX',p,flags);
cmd = sprintf('%s rmesh.cpp ../src/RectMeshUnstruct.cpp ../src/CodeAnalysis.cpp MexUtil.cpp',mc);
fprintf('%s\n', cmd);
eval(cmd);

% Configure options
isunix = 1;
use_omp = 1;

% Not sure I have all the flags right. The UNIX part works on my Ubuntu
% 64-bit 4-CPU system with use_omp = 0 or 1. At least one 32-bit Windows user
% has built the code successfully using this script.
p = '.';
flags = '-O -DHMMVP_MEXSVD';
if(use_omp) flags = [flags ' -DHMMVP_OMP ']; end
if(isunix)
  if(use_omp)
    flags = [flags [' CXXFLAGS="\$CXXFLAGS -fopenmp" ',...
		    'LDFLAGS="\$LDFLAGS -fopenmp"']];
  end
  flags = [flags ' CXXLIBS="\$CXXLIBS -Wall -lmwblas -lmwlapack"'];
  %mc = sprintf('mex -f %s/mexopts.sh -outdir %s -largeArrayDims %s',p,p,flags);
  mc = sprintf('mex -outdir %s -largeArrayDims %s',p,flags);
else
  % Definitely not sure about the Windows stuff.
  use_omp = 0;
  flags = [flags ' CXXLIBS="$CXXLIBS -lmwblas -lmwlapack'];
  if(~isempty(dir(sprintf('%s/mexopts.bat',p))))
    mc = sprintf('mex -f %s/mexopts.bat -outdir %s %s',p,p,flags);
  else
    mc = sprintf('mex -outdir %s %s',p,flags);
  end
end

mc = sprintf('%s -I../include -I../src -I../util/include -I../util/matlab -I.', mc);
eval(sprintf('%s hmcc_hd_mex.cpp ../src/Hd.cpp ../util/matlab/MexUtil.cpp',mc));
eval(sprintf('%s hm_mvp.cpp ../src/Hmat.cpp ../src/HmatIo.cpp',mc));

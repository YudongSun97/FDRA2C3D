%The file make.m didn't work for me. Instead, this command worked:

mex dc3dm_mex.cpp ../src/RectMeshUnstruct.cpp ../src/Triangulation.cpp ../src/MeshAnalyzer.cpp ../src/CodeAnalysis.cpp MexUtil.cpp -largeArrayDims -I.. -lmwblas -lmwlapack

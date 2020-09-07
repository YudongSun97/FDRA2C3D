HMMVP: Software to Compute Hierarchical-Matrix Matrix-Vector Products
Version 1.0
Andrew M. Bradley
ambrad@cs.stanford.edu
CDFM, Geophysics, Stanford
https://pangea.stanford.edu/research/CDFM/software/index.html


This package provides:
1. Pure Matlab routines to create an H-matrix approximation to a matrix;
2. C++ routines to compute matrix-vector products (MVP) and related operations
   using this H-matrix approximation.
3. A Matlab mex interface to these routines.

HMMVP is licensed as follows:
Open Source Initiative OSI - Eclipse Public License 1.0
http://www.opensource.org/licenses/eclipse-1.0


Matlab Installation
-------------------
You need a C++ compiler.
  If you're on a Unix-like or Mac system, install GNU g++.
  If you're on Windows, install either the cygwin or mingw ports of C++ (google
for "mingw"). Then obtain gnumex.m (google for "gnumex") and create a
mexopts.bat file pointing to the C++ compiler. Place the mexopts.bat file in
this directory. This page
    http://gnumex.sourceforge.net/#fortran_compile
provides a particularly clear description of how to get mex compilation working
on Windows for C, C++, and Fortran code.

Once Matlab is using the C++ compiler successfully,
1. At the Matlab command line, cd to the hmmvp directory.
2. Edit 'make.m'.
   a. If you are on a Unix-like or Mac system, set isunix = 1; if Windows, 0.
   b. If your computer has multiple CPUs, set use_omp = 1, otherwise 0.
3. Type 'make' at the command line. If you get no errors, you're done.
4. If 'make' fails and you had set use_omp = 1, set it to 0 instead and try
   again.
5. If 'make' still fails, email me with the error output.


Usage
-----
See 'ex.m' for examples of usage in Matlab.

The basic steps to compress the matrix B of Green's functions are as follows.
   1. Spatially decompose the mesh by calling hmcc_hd in Matlab. This is fast.
   2. Call hm('EstBfroLowBndCols',...) to estimate ||B||_F.
   3. Form the H-matrix by calling hm('Compress',...) in Matlab using the
      estimate from step 2. This can be time consuming for large matrices. Form
      it with a tolerance that is ~10 times tighter than what you really want.
   4. Form another H-matrix at a looser tolerance using the H-matrix from step
      3. This step produces an H-matrix that is just slightly more accurate than
      what is requested, which is ideal. Forming this second H-matrix is much
      faster than forming the first because the data from the first H-matrix is
      used extensively.
   5. Compute many matrix-vector products in Matlab using the mex wrapper hm_mvp
      or by calling the C++ code directly. While steps 1-4 must be done in
      Matlab and can be slow, this step is where the main computations of
      interest are performed. These can be performed in Matlab or in
      C++. Both OpenMP (shared memory) and MPI (distributed computers)
      implementations are available.

Matlab usage information is as follows:
  Compress: In Matlab, type 'help hmcc_hd' and 'help hm'.
  MVP: In Matlab, type 'help hm_mvp'.

The C++ interfaces to the spatial decomposition and MVP are complete, but
compression must be done in Matlab.
  OpenMP is used for all Hmat::Mvp() variants, but parallelization scales well
only for the full-matrix MVP at present: ie, B*x scales well, but B(rs,:)*x and
B(rs,cs)*x(cs) do not.
  MpiHmat.?pp implement parallelization of the full-matrix MVP in MPI. See
the header notes in MpiHmat.hpp for basic usage.
  The directory "messy" contains some C++ programs that illustrate C++ usage and
compilation. These are not meant to be compilable or runnable; but they may be
helpful in creating and building your own programs.

A very incomplete Fortran interface to just the full-matrix OpenMP-parallelized
MVP is available. The inteface is implemented in the files CFHmat.h and
CFHmat.cpp, and an example of using and building a Fortran 90+ file is
ex_mvp.f90.

You may find my report "H-Matrix and Block Error Tolerances"
    http://arxiv.org/abs/1110.2807
helpful when setting the error tolerance. In particular, you should specify an
error tolerance epsilon such that
    -log10(epsilon) - log10(cond(B)) >= (desired relative error in the MVP),
where B is the matrix. For example, if cond(B) = 1e5 and you want about 3 digits
of accuracy, set epsilon = 1e-8.

If you use this software in research that you publish, I would appreciate your
citing this software package as follows:

    A. M. Bradley, "HMMVP: Software to Form an H-Matrix and Compute
    H-Matrix-Vector Products", 2012.

Email me (ambrad@cs.stanford.edu) if you want to receive an email when I release
a new version.


Important version changes
-------------------------
Version 0.6:
1. Can't read previous versions' H-matrix files.
2. The interface to hmcc_hd has changed. It now takes in point clouds and
   returns permutations corresponding to the domain and range spaces. If the two
   spaces are the same, then pass the same point cloud for both arguments; and
   the two output permutations will be the same.
Version 0.8:
1. Can't read previous versions' H-matrix files.
Version 0.11:
1. Slight change to hm('Compress') call sequence.
Version 0.14:
1. What used to be Hmat<Real> is now TypedHmat<Real>. Hmat is now a superclass
   that hides the detail of whether the underlying matrix uses float or double
   type. There's a corresponding factory function, NewHmat, that lets the caller
   write
       Hmat* hm = NewHmat(filename, ...);
   The caller then doesn't have to know or specify whether file 'filename' uses
   float or double type. If to date you've used only the Matlab wrapper hm_mvp,
   this change does not affect your code.
Vertsion 0.15:
1. Subset MVP interface has changed.


Release log
-----------
0.0.  AMB. 29 Dec 2010. Early release.
0.1.  AMB. Compilation tweaks for Windows. (Thanks to Kaj Johnson.)
0.5.  AMB. 28 Jun 2011. First release.
0.6.  AMB. 27 Sep 2011. New release. Many changes.
0.7.  AMB. 3  Nov 2011. Hmat::FullBlocksIJS and hm_mvp('bblocks').
0.8.  AMB. 8  Nov 2011. Use an old H-matrix file in hm('Compress') to speed up
      creating the new one by using the option 'old_hmat_fn'.
0.9.  AMB. 16 Nov 2011. Hmat::SetThreads. Can now switch between number of
      threads dynamically. Also some small speedups.
0.10. AMB. 24 Nov 2011.
      - Fixed a bug that arose in pure C++ calls to Hmat<float>::Mvp when x and
        y are double.
      - Different OpenMP parallelization strategy in some routines for better
        speedup (wall-clock time / number of threads).
0.11. AMB. 12 Dec 2011.
      - MpiHmat.?pp. Full-matrix MVP in MPI.
      - Added directory "messy". This will contain non-buildable examples that
        may be of use.
0.12. AMB. 15 Dec 2011.
      - Added a very basic Fortran interface.
      - Improvements to MpiHmat.
0.13. AMB. 22 Mar 2012.
      - Fixed a bug in hm('FroNorm'), though the routine is not actually used.
      - Added a user-input check to hmcc_hd.m to detect duplicate points. (Thanks
        to Mark McClure.)
      - Added omp_set_num_threads(nthreads) to the start of all public
        entrypoints. This reestablishes the correct #threads in case another
	OpenMP-enabled library has changed it.
0.14. AMB. 18 May 2012.
      - TypedHmat<Real> is a subclass of Hmat. Now a container, such as
        vector<Hmat*>, can hold H-matrices having different (float or double)
	types.
      - As a consequence, hm_mvp no longer needs to be compiled for specifically
        float or double H-matrices.
      - Added the following public routines to hm.m:
        . BlackBoxCompress: Call this to get the best default options and to
          automatically 
        . BlocksPrecon: Get a decent preconditioner for the H-matrix.
        . SolveBxb: Use GMRES to solve B x = b.
        . CondEst1: Estimate ||B|_1.
	. NormInvEstFro: Estimate ||inv(B)||_F.
        . CondEstFro: Estimate ||B||_F.
        . New hm('Compress') option: Set 'tol_method' either to 'Bfro' (old
	  method) or 'BinvFro' (new method).
0.15. AMB.
      - Cleaned up some issues that made zero blocks not compress. Thanks Mark.
      - Added tol_method = 'abs'. This option lets the user specify a tolerance
        in whatever way is appropriate. The options 'Bfro' and 'BfroInv' are
        specializations addressing common cases.
      - hmcc_hd now returns a 4xNblocks double array rather than a struct
        array. I discovered that Matlab structs have too much extraneous memory
        associated with them for this application.
      - Added y = B(:,cs) x(cs).
      - Experimental: Mark McClure suggested the element-to-block map
        as an alternative implementation of the subset operations
           y(rs) = B(rs,cs) x(cs)
           y(rs) = B(rs,:)  x(:)
           y     = B(:,cs)  x(cs).
        The public interface is
           void TurnOnElementToBlockMap();
           void TurnOffElementToBlockMap();
        The map can add ~50% memory overhead and so cannot always be
        used. However, once the map is constructed, the algorithm is much faster
        than the default algorithm for very small subsets: those on the order of
        10s of columns or rows or smaller.
0.16. AMB.
      - Changed one part of the element-to-map block based on an idea by Mark
        McClure.
      - Changed how rs/cs_idxs are filled in when the E2B map is used.
0.17. AMB.
      - Namespace Hm -> hmmvp.
      - Everything is in hmmvp namespace now to avoid namespace pollution.
      - Hmat.hpp has been purified of the need to #include "Hmat.cpp" for
        template reasons.

To do
-----
- Implement efficiently scaling OpenMP and MPI implementations of B(rs,:)*x and
  B(rs,cs)*x(cs).
- Implement transpose subset products.
- MPI/C++ implementation of hm.m routines.


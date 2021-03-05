/* fdra2c 0.0: Fault Dynamics with a Radiation Damping Approximation
   MPI-enabled version
   AMB ambrad@cs.stanford.edu
   CDFM, Geophysics, Stanford University
   
   fdra2c is licensed as follows:
     Open Source Initiative OSI - Eclipse Public License 1.0
     http://www.opensource.org/licenses/eclipse-1.0

   make mode=p opt=-O
     parallelization by MPI
   make mode=s opt=-O
     serial with no MPI dependencies

   Parameters
   ----------
     nelem ncomp

     evolution = 'aging' or 'slip'
     mu0 a b d_c v0(s)
     [use_vcutoff vc_v1 vc_v2]

     eta(s) s_normal
     stress_fn =
       'ss': ss_k v_creep(s)
       'h-matrix': hm_filename [hm_scale (1)] hm_bc hm_symmetric [v_creep (1)]
         hm_symmetric = 1: hm_perm_p, hm_perm_q1, hm_perm_q2

     ti tf

     v_init chi_init [slip_init] [dlt_init]

     [stop_indicator] [stop_check_frequency]

     disp_every
     allow_overwrite
     [save_filename save_v_every save_slip_every save_state_every]
     [lineint_save_filename lineint_save_every lineint_wts]

   Notes
   -----
     hm_filename is the base name of the H-matrix files. These include:
       _comp11.hmat
       _comp12.hmat   (include these
       _comp21.hmat    ...
       _comp22.hmat    only if ncomp == 2)
       _compns.hmat   (optional normal-stress coupling component)

     If hm_hmat_symmetric, then permutation matrices (encoded as vectors q1,q2
     and p) Q and P must be supplied to map the domain to the H-matrix space (Q)
     and back (P). Indexing starts at 1.
       q has two columns since the symmetry requires two MVP. It maps the domain
     ordering to the H-matrix ordering: The H-matrix acts on xh in
         xh(i,j) = x(q(i,j)), j = 1,2.
       P injects the two columns y1 and y2 into the domain:
         y(i) = yh(p(i)), where yh = [y1; y2], y1 = B*xh(:,2), y2 = B*xh(:,2).
       For example, suppose an 8-element rectangular domain looks like this:
         [1 5
          2 6
          ---  <-- mirror symmetry across this line  
          3 7
          4 8]
     Vectorized, it becomes (1:8)'.
       The 4x8 H-matrix B was constructed so that it is applied as follows:
         y1 = B [2 1 3 4 6 5 7 8]'
         y2 = B [3 4 2 1 7 8 6 5]'
       To encode these transformations, set
         q(:,1) = [2 1 3 4 6 5 7 8]'
         q(:,2) = [3 4 2 1 7 8 6 5]'.
     The resulting full vector is assembed as
         vec([flipud(reshape(y1,2,2))
              reshape(y2,2,2)]).
     Hence the permutation applied to yh = [y1(:); y2(:)] to get
         y = yh(p)
     is
         p = [2 1 5 6 4 3 7 8].
     Because of the stacking in yh, 1:4 refer to the elements of y1; 5:8, of y2.

     hm_scale is a scalar multiple applied to both the H-matrix and the boundary
     condition. It's for convenience. It can account for changes in sign
     convention, problem scale, units, etc.

     The lineint_* options assume an nz x nx domain with fast index in the z
     direction and integration in the x direction. nz and nx are derived from
     the length of lineint_wts, which has length nx.  */

#include <stdio.h>
#include "util/include/Mpi.hpp"
#include "util/include/KeyValueFile.hpp"
#include "Fdra.hpp"
#include "HmatrixStressFn.hpp"
#include "MyCodeAnalysis.hpp"
using namespace fdra;
using namespace util;

static void Finalize (int code) {
  mpi::Finalize();
  exit(code);
}

int main (int argc, char** argv) {
  mpi::Init(argc, argv);
  bool am_root = mpi::AmRoot();

  if (argc != 2) {
    fprintf(stderr, "%s key-value-file\n", argv[0]);
    Finalize(-1);
  }

  KeyValueFile* kvf = NULL;
  if (am_root) {
    kvf = NewKeyValueFile();
    if (!kvf->Read(argv[1])) {
      fprintf(stderr, "Failed to read the key-value file %s.\n", argv[1]);
      DeleteKeyValueFile(kvf);
      Finalize(-1);
    }
  }
  Model* model = BuildModelFromKeyValueFile(kvf);
  if (!model) {
    if (am_root)
      fprintf(stderr, "The key-value file %s does not create a full model.\n",
              argv[1]);
    Finalize(-1);
  }

  if (!model->GetStressFn()) {
    // No default stress function was requested, so let's go with the H-matrix
    // one.
    model->SetStressFn(NewHmatrixStressFn(model, kvf));
  }

  if (am_root) DeleteKeyValueFile(kvf);

  if (!model || !model->IsOk()) {
    if (am_root)
      fprintf(stderr, "The key-value file %s does not create a full model.\n",
              argv[1]);
    DeleteModel(model);
    Finalize(-1);
  }

  if (!model->InitListeners()) {
    DeleteModel(model);
    Finalize(-1);
  }

  if (am_root) kvf = NewKeyValueFile();
  model->ToKvf(kvf);
  if (am_root) {
    kvf->Write("fdra.kvf");
    DeleteKeyValueFile(kvf);
  }

  mpi::Barrier(); Ca::GetTimer()->Tic(0);
  Go(model);
  mpi::Barrier(); Ca::GetTimer()->Toc(0);

  DeleteModel(model);

  caprint("pid %d: call et = %1.2e %1.2e factor = %1.2e | mvp %1.3e\n",
          mpi::Pid(), Ca::GetTimer()->TotEt(1), Ca::GetTimer()->TotEt(0),
          Ca::GetTimer()->TotEt(1)/Ca::GetTimer()->TotEt(0),
          Ca::GetTimer()->TotEt(11));

  Finalize(0);
  return 0;
}

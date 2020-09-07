function varargout = wrk(varargin)
  [varargout{1:nargout}] = feval(varargin{:});
end

function [ma me] = test_Ode()
  c.A = [1 1; -2 -1];
  if (true)
    c.A = kron(c.A, eye(15));
  end
  m = size(c.A,1);
  c.y0 = randn(m,1);
  c.ti = 0;
  c.tf = 10;
  c.reltol = 1e-4;
  c.abstol = 1e-6;
  c.initial_step = 0.1;
  
  fdra2c('WriteKvf', 'test_ode.ser', c);
  
  opts = odeset('reltol', c.reltol, 'abstol', c.abstol,...
                'initialstep', c.initial_step);
  t = tic();
  [ma.t ma.y] = ode23(@(t,y)c.A*y, [c.ti c.tf], c.y0, opts);
  etma = toc(t);

  function [yd serr] = OdeFn(t,dt,y)
    serr = false;
    yd = c.A*y;
  end
  function stop = OdeOutputFn(t,dt,y,msg)
    go.ssd = SaveStreamData('Write', go.ssd, [t(1); y(:)]);
    stop = false;
  end
  global go;
  go.ssd = SaveStreamData('Init', 'tmp.dat', 'nr', m + 1);
  opts.OutputFcn = @OdeOutputFn;
  t = tic();
  ode('Pair23', @OdeFn, [c.ti c.tf], c.y0, opts);
  etmem = toc(t);
  SaveStreamData('Finalize', go.ssd);
  d = SaveStreamData('Read', 'tmp.dat');
  mem.t = d(1,:);
  mem.y = d(2:end,:);
  
  np = min(m, 8);
  fprintf(1, 'np = %d\n', np);
  t = tic();
  [s r] = system(sprintf('mpirun -np %d ./test_Ode;', np));
  etme = toc(t);
  d = SaveStreamData('Read', 'test_ode.dat');
  me.t = d(1,:);
  me.y = d(2:end,:);
  
  fprintf(1, 'et mat, mem, me: %1.3e %1.3e %1.3e\n', etma, etmem, etme);
  
  sp(311); plot(ma.t, ma.y, 'b.-', me.t, me.y, 'ro-', mem.t, mem.y, 'gx-');
  sp(312); plot(ma.t(1:end-1), diff(ma.t), 'b.-',...
                me.t(1:end-1), diff(me.t), 'r.-',...
                mem.t(1:end-1), diff(mem.t), 'gx-');
  sp(313); sy(mem.t(1:end-1), abs(diff(me.t(:)) - diff(mem.t(:))), '.-');
end

function [c s n] = test_ss(nrep)
  cfn = '/scratch/ambrad/fdra2c/mcode';
  c = AssembleStruct('fs_ss', cfn);
  
  mt = tic();
  fdra('fs_ss');
  em = toc(mt);
  fprintf(1, 'mat et = %1.2e\n', em);
  s = Util('qloadr', c.saveFn, 2);

  [n et] = test_ss_cpp(c, nrep);
  fprintf('speedup = %1.2f\n', em/et);
  
  clf;
  plot(s.psi,s.gam,'bo',...
       log(n.v/s.c.v0(1)),log(n.theta*(s.c.v0(1)/s.c.d_c(1))),'r.');
end

function c = AssembleStruct(s, saveFn)
  try
    % Save the script so we can refer to it later
    if (~strcmp(s(end-1:end),'.m'))
      s = [s '.m'];
    end
    fid = fopen(s,'r');
    script = fread(fid,inf,'char=>char')';
    fclose(fid);
    % Run the script
    eval(script);
  catch
    lasterr
    script = [];
    eval(s);
  end

  save foo.mat;
  c = load('foo');
end

function bds = Segment(n, nthreads)
  bds = zeros(1, nthreads + 1);
  k = floor(n / nthreads);
  extra = n - k*nthreads;
  b = 0;
  for (i = 1:nthreads)
    bds(i) = b;
    b = b + k;
    if (extra > 0)
      b = b + 1;
      extra = extra - 1;
    end
  end
  bds(nthreads+1) = n;
end

function WriteSsInput(s,nrep)
% s is saved variables from a fdra script.
% nrep makes nrep replicas of the problem to test parallelism.
  c.ncomp = 1;
  if (numel(s.a) == 1)
    o = ones(1, nrep);
    c.nelem = nrep;
  else
    o = 1;
    c.nelem = numel(s.a);
  end

  if (s.evolution == 1) c.evolution = 'aging';
  else c.evolution = 'slip'; end
  c.mu0 = s.mu_0;
  c.b = s.b*o;
  c.a = s.a*o;
  c.d_c = s.Dc*o;
  c.v0 = s.V0;
  c.use_vcutoff = s.use_vcutoff;
  c.vc_v1 = s.v1*o;
  c.vc_v2 = s.v2*o;

  year2sec = s.year2sec;
  sec2year = 1/year2sec; 
  c.v_creep = s.Vpl_creep;
  c.stress_fn = 'ss';
  c.eta = s.eta;
  c.s_normal = s.s_normal*o;
  c.ss_k = s.ss_k*o;

  c.ti = s.ts;
  c.tf = s.tend;
  
  c.v_init = s.v_init*o;
  c.chi_init = s.chi_init*o;
  
  c.stop_indicator = 'stop.ind';
  c.stop_check_frequency = 1;
  c.disp_every = 1000;
  
  c.save_filename = '/scratch/ambrad/fdra2c/ss';
  c.allow_overwrite = 1;
  se = s.Nsavepoints;
  c.save_v_every = se;
  c.save_slip_every = se;
  c.save_state_every = se;
  
  fdra2c('WriteKvf', '/scratch/ambrad/fdra2c/ss.kvf', c);
end

function [s et] = test_ss_cpp(c,nrep)
  WriteSsInput(c,nrep);
  if (numel(c.a) > 1) nrep = numel(c.a); end
  t = tic();
  [~,r] = system(sprintf('mpirun -np %d ./fdra /scratch/ambrad/fdra2c/ss.kvf',...
                         min(4, nrep)));
  et = toc(t);
  fprintf(1, 'mpi et = %1.2e\n', et);
  r
  bfn = '/scratch/ambrad/fdra2c/ss';
  v = SaveStreamData('Read', [bfn '_v.sdf'], [], 1:3);
  slip = SaveStreamData('Read', [bfn '_slip.sdf'], [], 1:3);
  theta = SaveStreamData('Read', [bfn '_theta.sdf'], [], 1:3);
  clf;
  if (true)
    sp(311); plot(v(1,:), log10(v(2:end,:)), '.-'); ylabel('log10 v');
    sp(312); plot(slip(1,:), slip(2:end,:), '.-'); ylabel('slip');
    sp(313); plot(theta(1,:), log10(theta(2:end,:)), '.-'); ylabel('log10 \theta');
  else
    plot(log10(v(2:end,:)), log10(theta(2:end,:)), '.-');
  end
  if (relerr(v(2,:), v(end,:)) > 0) kb; end
  s.t = v(1,:);
  s.v = v(2,:);
  s.slip = slip(2,:);
  s.theta = theta(2,:);
end

function test_physics()
 v0  = 1.00000e-06;
 mu0 = 6.00000e-01;
 a   = 1.62000e-02;
 b   = 1.80000e-02;
 d_c = 1.00000e-04;

 t     = 2.00000e+00;
 slip  = 1.00000e-04;
 v     = 3.00000e-01;
 theta = 5.00000e-02;
 
 k = 1.8e7;
 v_creep = 1.267523512561158e-09;
 tau = k*(v_creep*t - slip)
 
 c.v0 = v0;
 c.mu_0 = mu0;
 c.a = a;
 c.b = b;
 c.d_c = d_c;
 c.evolution = 2;
 c.use_flash = 0;
 c.use_vcutoff = 0;
 
 psi = log(v/v0);
 gamma = log(theta*v0/d_c);
 [mu mu_psi mu_gamma] = fdra2('Friction', c, psi, gamma)
 
 thd_o_th_el2 = fdra2('Evolve', c, psi, gamma)
 c.evolution = 1;
 thd_o_th_el1 = fdra2('Evolve', c, psi, gamma)
end

function c = test_hm_sym(savefn, ow)
  switch (3)
   case 1
    tol = 3;
    pfn = 'pnx88xnz121';
   case 2
    tol = 2;
    pfn = 'pnx670nz2117';
   case 3
    tol = 2;
    pfn = 'pnx1337nz2478';
  end
  load(['evolve/' pfn '.mat']);
  p.d_c = p.d_c;
  p = evolve_test('Reflect', p);
  flds = {'nx' 'nz' 'z' 'zc' 'x' 'xc'};
  for (i = 1:length(flds)) c.(flds{i}) = p.(flds{i}); end
  c.nelem = numel(p.a);
  c.ncomp = 1;
  c.v0 = p.init.v0;
  c.v_creep = p.init.v_creep;
  c.use_flash = 0;
  c.use_vcutoff = 1;
  c.evolution = 'aging';
  c.eta = 1e-6;
  
  c.hm_filename = sprintf('/scratch/ambrad/fdra2c/bem_id1tol-%d.0nx%dnz%dsym',...
                          tol, p.nx, p.nz/2);
  c.hm_symmetric = 1;
  c.hm_scale = 1;
  
  flds = {'a' 'b' 'd_c' 's_normal'};
  for (i = 1:numel(flds)) c.(flds{i}) = p.(flds{i})(:); end
  c.vc_v1 = p.v1(:);
  c.vc_v2 = p.v2(:);
  v_match = 1e-12;
  if (c.use_vcutoff)
    c.mu0 = 0.6 + (c.a - c.b).*log(v_match/c.v0) +...
            c.a.*log(c.vc_v1/v_match + 1) - c.b.*log(c.vc_v2/v_match + 1);
  else
    c.mu0 = 0.6*ones(size(c.a));
  end

  lo = load([c.hm_filename '_p.mat']);
  n = size(p.s_normal, 1)/2;
  i = n + find(p.s_normal(n+1:end,1) > 2*p.s_normal(n+1,1), 1);
  z_tran = p.z(i+1);
  if (true)
    symm = 1;
    switch (3)
     case 1 % bc1
      c.hm_bc = evolve_test('SetupBCs', lo.p, z_tran);
      c.hm_bc = c.hm_bc(1:length(c.hm_bc)/2);
     case 2 % bc2
      c.hm_bc = fdra2i('SetupBCs', lo.p, 1, 1, 1, 0);
      c.hm_bc = c.hm_bc(:,1) + c.hm_bc(:,2) + c.hm_bc(:,3);
     case 3 % bc3
      symm = 0;
      xlim = 121e3;
      p.mu = lo.p.mu;
      p.nu = lo.p.nu;
      c.hm_bc = evolve_test('SetupAllSidesBCs', p, xlim);
    end
    if (symm)
      c.hm_bc = reshape(c.hm_bc, p.nz/2, p.nx);
      c.hm_bc = [c.hm_bc(end:-1:1,:); c.hm_bc];
      c.hm_bc = c.hm_bc(:);
    end
    bc = c.hm_bc;
    save foo bc;
  else
    pr('loading bc from file\n');
    load foo;
    c.hm_bc = bc;
  end

  vec = @(x) x(:);
  nzh = p.nz/2;
  is = reshape(1:p.nz*p.nx, p.nz, p.nx);
  is = [vec(is(nzh:-1:1,:)); vec(is(nzh+1:end,:))];
  c.hm_perm_q1 = is;
  nh = numel(is)/2;
  c.hm_perm_q2 = is([nh+1:end 1:nh]);
  is = reshape(1:nzh*p.nx, nzh, p.nx);
  is = [flipud(is); nh + is];
  c.hm_perm_p = is(:);
  
  [psi gam chi] = evolve_test('InterpICs',...
      p.xc*1e-3,'~/dd/evolve/ffref_dc15_ics_89281');
  v = p.init.v0*exp(psi);
  c.v_init = vec(repmat(v(:)',size(p.a,1),1));
  gamma_init = vec(repmat(gam(:)',size(p.a,1),1));
  psi = log(c.v_init/p.init.v0);
  % Mess with N, S v.s. areas.
  switch (1)
   case 1
    % Set the boundary area to plate rate. This seems to set off
    % accelerations in the v.s. area at the bdies in the transition region,
    % but it may be that they stay controlled and then go away.
    m = p.a > p.b;
    psi(m) = log(c.v_creep/p.init.v0);
   case 2
    % Make the boundary slip speed fall off away from the E boundary. Bring
    % it down to 1e-3 plate rate by the transition region. This is causing
    % spiky behavior in the N, S bdy regions.
    x_cutoff = 1.199e5;
    zlim = p.zc(find(p.zc(:) > 0 & p.a(:,1) >= p.b(:,1),1));
    m = p.a > p.b & abs(p.Zc) > zlim;
    psi = reshape(psi, size(m));
    o = ones(c.nz,1);
    a = p.xc - x_cutoff;
    a = a/a(end);
    a(a < 0) = 0;
    a = 1e-3*(1 - a) + 1*a;
    q = o*log(a*c.v_creep/c.v0);
    % Choose the smaller of the IC curve and this artificial IC.
    psi(m) = min(psi(m), q(m));
   otherwise
    % Don't modify the boundary, but don't let anything be above plate rate.
    psi(c.v_init > c.v_creep) = log(c.v_creep/c.v0);
  end
  % Add roughness.
  m = p.Xc >= 130e3 & p.a < p.b;
  psi(m) = psi(m) + 1*randn(size(psi(m)));
  c.slip_init = zeros(size(c.v_init));
  gamma_init(m) = -psi(m);
  c.chi_init = psi(:) + gamma_init(:);
  c.v_init = p.init.v0*exp(psi(:));
  
  if (true)
    % Radically reshape s_normal to try to make v as non-spiky as possible in
    % the boundary region. [Later:] Nope, not really helping. I think I just
    % need to make the resolution better.
    c.s_normal = repmat(p.s_normal(end/2,:), c.nz, 1);
    c.s_normal = c.s_normal(:);
    w = c.a - c.b;
    w(w<0) = 0;
    w = w/max(w(:));
    m = c.a > c.b;
    c.s_normal(m) = (1-w(m)).*c.s_normal(m) + w(m)*max(c.s_normal(:));
  end
    
  c.ti = 0;
  c.tf = 1e12;
  
  c.disp_every = 5;
  c.stop_indicator = 'stop.ind';
  c.stop_check_frequency = 5;
  c.save_filename = savefn;
  c.save_v_every = 100;
  c.save_slip_every = 500;
  c.save_state_every = 500;
  c.lineint_save_filename = c.save_filename;
  c.lineint_save_every = 5;
  is = find(p.a(end/2,:) < p.b(end/2,:));
  is(end+1) = is(end) + 1;
  c.lineint_wts = zeros(1, c.nx);
  c.lineint_wts(is(1:end-1)) = diff(p.x(is))/(p.x(is(end)) - p.x(is(1)));
  c.allow_overwrite = ow;

  fdra2c('WriteKvf', [savefn '.kvf'], c, c.allow_overwrite);
end

function t = qload_t(fn, stride)
  t = SaveStreamData('Read', [fn '_v.sdf'], stride, 1, 1);
end

function s = qload_sym(fn, stride, do_stride, varargin)
  [want_slip want_c] = process_options(...
      varargin, 'want_slip', false, 'want_c', false);
  if (nargin < 3) do_stride = numel(stride) == 1; end
  c = fdra2c('ReadKvf', [fn '.kvf']);
  flds = {'nx' 'nz' 'x' 'z' 'xc' 'zc'};
  for (i = 1:length(flds)) s.(flds{i}) = c.(flds{i}); end
  s.v = SaveStreamData('Read', [c.save_filename '_v.sdf'],...
                       stride, [], do_stride);
  s.t = s.v(1,:);
  s.v = s.v(2:end,:);
  if (want_slip)
    s.slip = SaveStreamData('Read', [c.save_filename '_slip.sdf'],...
                            stride, [], do_stride);
    s.slip = s.slip(2:end,:);
    s.theta = SaveStreamData('Read', [c.save_filename '_theta.sdf'],...
                             stride, [], do_stride);;
    s.theta = s.theta(2:end,:);
  end
  s.fh = @(I) reshape(I, s.nz, s.nx);
  s.fhi = @(I) fdra2i('Interp', s, I);
  if (want_c) s.c = c; end
end

function s = qload_sym_lix(fn, stride, do_stride)
  if (nargin < 3) do_stride = numel(stride) == 1; end
  c = fdra2c('ReadKvf', [fn '.kvf']);
  s = [];
  try
    s.v = SaveStreamData('Read', [c.lineint_save_filename '_lix.sdf'],...
                         stride, [], do_stride);
    s.t = s.v(1,:);
    s.v = s.v(2:end,:);
    s.zc = c.zc;
  catch
    error('Failed to read lineint file.');
  end
end

function s = interp_lix(s, fac)
  if (nargin < 2) fac = 1; end
  t = linspace(s.t(1), s.t(end), fac*numel(s.t));
  s.v = interp1(s.t, s.v', t)';
  s.t = t;
  zc = linspace(s.zc(1), s.zc(end), 1.5*fac*numel(s.zc));
  s.v = interp1(s.zc, s.v, zc);
  s.zc = zc;
  s.iml = @() imagesc(s.t/(3600*24*365.25),s.zc*1e-3,log10(s.v));
end

function test_hm_vanilla()
  if (true)
    nx = 100;
    ny = 25;
    pad = 7;
  else
    nx = 400;
    ny = 125;
    pad = 28;
  end
  [I J] = meshgrid(1:nx, 1:ny);
  
  c.nelem = ny*nx;
  c.ncomp = 1;
  c.v0 = 1e-6;
  c.v_creep = 1e-9;
  c.use_vcutoff = 1;
  c.evolution = 'aging';
  c.eta = 1e-6;
  
  o = ones(ny,nx);
  b = 0.018;
  c.b = b*o;
  c.a = 0.9*c.b;
  c.a(I <= pad | I > nx - pad) = 1.1*b;
  c.a(J <= pad | J > ny - pad) = 1.1*b;
  c.vc_v2 = 1e-6*o;
  c.vc_v1 = 1e2*o;
  c.s_normal = 100e6*o;
  c.d_c = 1e-4*o;
  v_match = 1e-12;
  c.mu0 = 0.6 + + (c.a - c.b).*log(v_match/c.v0) +...
          c.a.*log(c.vc_v1/v_match + 1) - c.b.*log(c.vc_v2/v_match + 1);
  
  c.chi_init = 0*o;
  c.v_init = c.v_creep*o;

  c.stress_fn = 'h_matrix';
  c.hm_filename = '/scratch/ambrad/fdra2c/Hmat_vanilla_small';
  c.hm_bc = o;
  c.hm_scale = 1;
  c.hm_symmetric = 0;
  
  c.allow_overwrite = 1;
  c.save_filename = '/scratch/ambrad/fdra2c/vanc';
  c.save_v_every = 50;
  c.save_slip_every = 200;
  c.save_state_every = 400;
  c.stop_check_frequency = 10;
  c.stop_indicator = 'stop.ind';
  c.disp_every = 1000;
  
  c.ti = 0;
  c.tf = 1e9;
  
  c.nx = nx;
  c.ny = ny;
  fdra2c('WriteKvf', '/scratch/ambrad/fdra2c/van.kvf', c);
end

function s = qload_vanilla(stride)
  s.c = fdra2c('ReadKvf', '/scratch/ambrad/fdra2c/van.kvf');
  s.v = SaveStreamData('Read', [s.c.save_filename '_v.sdf'], stride, [], 1);
  s.t = s.v(1,:);
  s.v = s.v(2:end,:);
  s.slip = SaveStreamData('Read', [s.c.save_filename '_slip.sdf'], stride, [], 1);;
  s.slip = s.slip(2:end,:);
  s.fh = @(I) reshape(I, s.c.ny, s.c.nx);
end

function call_fdra2i(c, nthreads)
% c is the kvf struct.
  if (nargin < 2) nthreads = 4; end
  basefn = c.hm_filename;
  savefn = c.save_filename;
  s.pure_ode = 1;
  s.use_nthreads = nthreads;
  s.monitor_threads = 0;
  s.disp_text = c.disp_every;
  if (strcmp(c.evolution, 'slip')) s.evolution = 2; else s.evolution = 1; end
  s.v1 = c.vc_v1(:);
  s.v2 = c.vc_v2(:);
  s.Chyd = [];
  c.mu0 = 0.6;
  s.v_init = c.v_init(:);
  s.gamma_init = c.chi_init(:) - log(c.v_init(:)./c.v0);
  flds = {'nelem' 'ncomp' 'v0' 'v_creep' 'use_vcutoff' 'mu0' 'eta'};
  for (i = 1:length(flds)) s.(flds{i}) = c.(flds{i}); end
  flds = {'a' 'b' 'd_c' 's_normal'};
  for (i = 1:length(flds)) p.(flds{i}) = c.(flds{i})(:); end
  if (isfield(c, 'ny'))
    p.nz = c.ny;
  else
    p.nz = c.nz;
  end
  p.nx = c.nx;
  ow = c.allow_overwrite;
  Bbc = c.hm_bc(:);
  if (c.hm_symmetric)
    symm_type = 2;
    Bbc = reshape(Bbc, p.nz, p.nx);
    Bbc = Bbc(p.nz/2+1:end,:);
    Bbc = Bbc(:);
  else
    symm_type = 0;
  end
  save([c.hm_filename '_bc.mat'], 'Bbc');
  fdra2i('Run', basefn, savefn, s, 'parms', p, 'symm_type', symm_type,...
         'overwrite', ow);
end

function test_perms()
  n = 10;
  p1 = randperm(n);
  p2 = randperm(n);
  is = 1:n;
  is = is(p1);
  is = is(p2);
  is21a = is;
  is(p2) = is;
  is(p1) = is;
  is
  p21 = p1(p2);
  is21b = is(p21);
  is(p21) = is21b;
  is
end

function d2 = dd_new(x, y)
  d2.x = x;
  if (nargin > 1) d2.y = y; else d2.y = 0; end
end

function c = dd_add(a, b)
  if (isstruct(b))
    % Adapted from ddadd in David H. Bailey's DDFUN. This implementation is
    % based on older work. In particular:
    t1 = a.x + b.x;
    if (true)
      %   Two-Sum(a.x, b.x) in Shewchuk's "Robust Geometric Predicates" paper,
      % attributed to one of D.E. Knuth's books. Two-Sum is a branch-less
      % alternative to a branching call to Fast-Two-Sum(a.x, a.b) based on |a.x|
      % >= |b.x|.
      e = t1 - a.x;
      t2 = ((b.x - e) + (a.x - (t1 - e)));
    else
      %   Branch-based Two-Sum.
      % Fast-Two-Sum( . , . )
      if ((a.x > b.x) == (a.x > -b.x))
        % Fast-Two-Sum(a.x , b.x)
        t2 = b.x - (t1 - a.x);
      else
        % Fast-Two-Sum(b.x , a.x)
        t2 = a.x - (t1 - b.x);
      end
    end
    % Now accumulate the low-order parts.
    t2 = t2 + a.y + b.y;
    %   Fast-Two-Sum(t1, t2) in Shewchuk's paper, attributed to
    %     T.J. Dekker 1971.
    % Applicable because |t1| >= |t2|.
    c.x = t1 + t2;
    c.y = t2 - (c.x - t1);
  else
    % Specialized to the case b.y = 0.
    t1 = a.x + b;
    e = t1 - a.x;
    t2 = ((b - e) + (a.x - (t1 - e))) + a.y;
    c.x = t1 + t2;
    c.y = t2 - (c.x - t1);
  end
end

function d2 = dd_mult(d2, d)
  if (~isstruct(d))
    d2.x = d*d2.x;
    d2.y = d*d2.y;
  else
    error('not impl');
  end
end

function d1 = dd_tod1(d2)
  d1 = d2.x;
end

function test_ddouble(testno)
  digits 64;
  switch (testno)
   case 1
    t1 = pi;
    t2 = exp(1)*1e-12;
    t3 = -pi;
    t4 = 1e-25;
    t5 = 20;
    a = ((((vpa(t1) + vpa(t2)) + vpa(t3)) + vpa(t4)) + vpa(t5)) - vpa(t5)
    d2 = dd_new(t1); d2
    d2 = dd_add(d2, t2); d2
    d2 = dd_add(d2, t3); d2
    d2 = dd_add(d2, t4); d2
    d2 = dd_add(d2, t5); d2
    d2 = dd_add(d2, -t5); d2
    d1 = dd_tod1(d2);
    av = double(a);
    d = ((((t1 + t2) + t3) + t4) + t5) - t5;
    pr('%1.15e %1.15e %1.15e %1.3e %1.3e\n', d1, av, d, relerr(av, d1),...
       relerr(av, d));
   case 2
    t1 = pi;
    t2 = exp(1)*1e-12;
    t3 = -pi;
    a = vpa(t1) + vpa(t2);
    am = a * -1.01;
    a = a + am;
    d2 = dd_new(t1); d2
    d2 = dd_add(d2, t2); d2
    d2m = dd_mult(d2, -1.01); d2m
    d2 = dd_add(d2, d2m); d2
    d1 = dd_tod1(d2m);
    amv = double(am);
    fprintf(1,'%1.15e %1.15e %1.3e\n', d1, amv, relerr(amv, d1));
   case 3
    t = dd_new(1.0);
    h = 1e-18;
    tnew = dd_add(t, h);
    h = dd_tod1(dd_add(tnew, dd_mult(t, -1.0)))
  end
end

function s = quick_view(fn, tstride)
  t = wrk('qload_t', fn, 1); len(t)
  s = wrk('qload_sym', fn, len(t), 0, 'want_c', 1);
  fi(1); clf; iml(s.fh(s.v(:,end))); cb; ca = caxis(); caxis([-12 ca(2)]);
  fi(2); quick_view_lix(fn, tstride, s.c.z);
  fi(3); clf; sp(7,1,1:4);
  [I iz ix] = s.fhi(s.v(:,end));
  imagesc(ix*1e-3, iz*1e-3, log10(I));
  ca = caxis(); caxis([-12 ca(2)]); axis equal; axis image;
  d_c = s.fh(s.c.d_c);
  sp(7,1,5); plot(s.xc*1e-3, d_c(end/2,:)); axis tight;
  s_normal = s.fh(s.c.s_normal);
  sp(7,1,6); plot(s.xc*1e-3, s_normal(end/2,:)*1e-6); axis tight;
  amb = s.fh(s.c.a - s.c.b);
  sp(7,1,7); plot(s.xc*1e-3, amb(end/2,:)); axis tight;
  s = rmfield(s, 'c');

  fi(4);
  t = SaveStreamData('Read',[fn '_slip.sdf'],1,1,1); len(t)
  slip = SaveStreamData('Read', [fn '_slip.sdf'], len(t), [], 0);
  slip = s.fh(slip(2:end));
  sp(211); imagesc(slip); cb;
  sp(212); plot(s.xc*1e-3,slip(end/2,:));
end

function si = quick_view_lix(fn, tstride, z)
  si = wrk('qload_sym_lix', fn, tstride, 1);
  zw = diff(z)/(z(end) - z(1));
  si.v = si.v.*repmat(zw(:), 1, size(si.v,2));
  s2y = 1/(365.25*24*3600);
  clf;
  sp(212); plot(si.t(2:end)*s2y, diff(si.t)/tstride, '.-');
  sp(211);
  try
    si = wrk('interp_lix', si, 2);
    si.iml(); cb;
  catch
    iml(si.v); cb;
  end
end

function q = PlotFractionHstar(c)
  mu = 3e10;
  nu = 0.25;
  mu_hat = mu/(1 - nu);
  hs = pi*mu_hat*c.d_c ./ (4*c.s_normal.*(c.b - c.a));
  hs(c.a >= c.b) = nan;
  [dZ dX] = meshgrid(diff(c.x), diff(c.z));
  mD = max(dZ, dX);
  q = reshape(hs, c.nz, c.nx)./mD;
  imagesc(q);
  cb; title('h^*/(max dim)'); caxis([0 20]);
end

function scratch()
  s=wrk('qload_sym_lix',fn,1,1);clf;sp(211);iml(s.v);cb;caxis([-10 -9.4]);sp(212);plot(s.t(2:end)*s2y,diff(s.t),'.-')
   t=wrk('qload_t',fn,1);len(t),s=wrk('qload_sym',fn,len(t),0);clf;iml(s.fh(s.v(:,end)));cb;ca=caxis();caxis([-12 ca(2)])

   t=SaveStreamData('Read',[fn '_slip.sdf'],1,1,1);len(t),slip=SaveStreamData('Read',[fn '_slip.sdf'],len(t),[],0);slip=s.fh(slip(2:end));sp(211);imagesc(slip);sp(212);plot(s.xc*1e-3,slip(end/2,:));
   
   [X Z]=meshgrid(p.x,p.z);
   imagesc(max(diff(Z(:,1:end-1)),diff(X(1:end-1,:)')'))
end

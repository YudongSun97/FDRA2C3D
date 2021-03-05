function varargout = wrk(varargin)
  [varargout{1:nargout}] = feval(varargin{:});
end

function c = WriteKvf(N, tol, prob, c_use, varargin)
  if (nargin < 4) c_use = []; end
  [nrepeat] = process_options(varargin, 'nrepeat', 10);
  switch (prob)
   case 1
    p = ex('Parameters', 'delta', 1e-4, 'N', N, 'tol', tol, 'wiggle', 0.0,...
           'geom', 'surfcube');
    close;
    c.greens_fn = 'inverse-r';
    c.X = p.X;
    pr('N = %d\n', size(p.X,2));
    c.order = 3;
    c.delta = 0;
   case 2
    c.greens_fn = 'okada-rect-tensor-mesh';
    m = round(sqrt(N/2));
    c.x = linspace(-3, 3, 2*m);
    c.eta = linspace(0, 1, m);
    if (true)
      % Add some noise to get away from innocuous, but annoying for testing,
      % instability in Hd for uniform meshes.
      dx = c.x(2) - c.x(1);
      c.x = c.x + 0.01*dx*randn(size(c.x));
      deta = c.eta(2) - c.eta(1);
      c.eta = c.eta + 0.01*deta*randn(size(c.eta));
    end
    pr('N = %d\n', 2*(m - 1)^2);
    c.dipdeg = -12; c.depth_min = 0; c.y_min = 0;
    c.mu = 3e10; c.nu = 0.25;
    c.disl_dip = 1; c.disl_strike = 0; c.disl_tensile = 0;
   case 3
    c.greens_fn = 'okada-subduct-symm-rect-periodic-tensor-mesh';
    m = round(sqrt(N/2));
    c.x = linspace(-3, 3, 2*m+1);
    c.x = c.x(m+1:end);
    c.eta = linspace(0, 1, m);
    if (false)
      dx = c.x(2) - c.x(1);
      c.x = c.x + 0.01*dx*randn(size(c.x));
      c.x(1) = 0;
      deta = c.eta(2) - c.eta(1);
      c.eta = c.eta + 0.01*deta*randn(size(c.eta));
    end
    pr('N = %d\n', 2*(m - 1)^2);
    c.dipdeg = -12; c.depth_min = 0; c.y_min = 0;
    c.mu = 3e10; c.nu = 0.25;
    c.nrepeat = nrepeat;
  end
  c.allow_overwrite = true;
  bfn = sprintf('/scratch/ambrad/test/test%d%d',...
                log10(tol), ~isempty(c_use));
  c.write_hmat_filename = [bfn '.hmat'];
  c.bc_filename = [bfn '.bc'];
  if (isempty(c_use))
    c.write_hd_filename = '/scratch/ambrad/test/test.hd';
  else
    if (isfield(c_use, 'write_hd_filename'))
      c.use_hd_filename = c_use.write_hd_filename;
    else
      c.use_hd_filename = c_use.use_hd_filename;
    end
    c.use_hmat_filename = c_use.write_hmat_filename;
  end
  c.tol = tol;
  kvf('Write', '/scratch/ambrad/test/test.kvf', c, true);
end

function basefn = RunFdra2iSetup(c, Bfro, symm, dip_only)
  p.ncomp = 1;
  p.x = c.eta - c.eta(1);
  p.z = c.x;
  p.nx = numel(p.x) - 1;
  p.nz = numel(p.z) - 1;
  p.dip = -c.dipdeg;
  basefn = fdra2i('Setup', '/scratch/ambrad/test/', p, 'tol', c.tol,...
                  'overwrite', 1, 'Bfro', Bfro, 'dip_only', dip_only,...
                  'symm', symm);
end

function o = CmpToHm(N, tol, ncpu, prob)
  estr = sprintf(['mpirun -np %d ./bin/hmmvpCompress ',...
                  '/scratch/ambrad/test/test.kvf'], ncpu);
  pr('%s\n', estr);
  c = WriteKvf(N, tol, prob); tic(); system(estr); o.et1 = toc();
  pr('et: %f\n', o.et1);
  in = hm('HmatInfo', c.write_hmat_filename);
  o.Bfro = in.tol / tol;
  o.c = c;
  switch (prob)
   case 2, symm = 0; dip_only = 0;
   case 3, symm = 1; dip_only = 1;
  end
  t = tic(); o.bfn = RunFdra2iSetup(c, o.Bfro, symm, dip_only); o.et2 = toc(t);
  pr('et: %f\n', o.et2);
  if (N <= 3000)
    o.Bn = Extract(c.write_hmat_filename);
    o.Bh = Extract([o.bfn '_comp11.dat']);
    sp(221); iml(o.Bn); cb;
    sp(222); iml(o.Bh); cb;
    sp(234); iml((o.Bn-o.Bh)./o.Bh); cb;
    pr('nwre = %e\n', relerr(o.Bh, o.Bn));
  end
  o.cfn = GetCompressFactor(c.write_hmat_filename);
  o.chn = GetCompressFactor([o.bfn '_comp11.dat']);
  sp(235); o.bcn = ShowBc(o.c);
  lo = load([o.bfn '_bc.mat']);
  o.bch = reshape(lo.Bbc(1:end/2), len(o.c.x)-1, len(o.c.eta)-1);
  sp(236); iml(o.bch); cb;
  pr('bc nwre = %e\n', relerr(o.bch, o.bcn));
  pr('%f vs %f for speedup %1.1f\n', o.et1, o.et2, o.et2/o.et1);
end

function wtf(o)
  function same = Same(b1, b2)
    same =...
        b1.r0 == b2.r0 && b1.c0 == b2.c0 && b1.m == b2.m && b1.n == b2.n &&...
        relerr(b1.B, b2.B) <= 10*eps &&...
        size(b1.U,2) == size(b2.U,2) && relerr(b1.U, b2.U) <= 10*eps &&...
        relerr(b1.V, b2.V) <= 10*eps;
  end
  
  [n.bs n.p n.q n.nnz n.m n.n] = hm('ReadBlocks', o.c.write_hmat_filename);
  [h.bs h.p h.q h.nnz h.m h.n] = hm('ReadBlocks', [o.bfn '_comp11.dat']);
  if (relerr(h.p, n.p) > 0 || relerr(h.q, n.q) > 0)
    'hrm', kb;
  end
  for (i = 1:numel(n.bs))
    if (~Same(n.bs(i), h.bs(i)))
      kb
    end
  end
end
  
function cf = GetCompressFactor(fn)
  in = hm('HmatInfo', fn);
  pr('realp = %d\n', in.realp);
  [id nnz] = hm_mvp('i', fn);
  hm_mvp('c', id);
  cf = in.m*in.n/nnz;
  pr('compress: %1.3f\n', cf);
end

function B = Extract(fn)
  in = hm('HmatInfo', fn);
  pr('realp = %d\n', in.realp);
  [id nnz] = hm_mvp('i', fn);
  pr('compress: %1.3f\n', in.m*in.n/nnz);
  B = hm_mvp('extract', id, 1:in.m, 1:in.n);
  hm_mvp('c', id);
end

function bc = ShowBc(c)
  fid = fopen(c.bc_filename, 'r');
  bc = fread(fid, inf, 'double');
  fclose(fid);
  bc = reshape(bc, numel(c.x) - 1, numel(c.eta) - 1);
  iml(bc); cb;
end

function c = test_compress(N, ncpu)  
  estr = sprintf(['mpirun -np %d ./bin/hmmvpCompress ',...
                  '/scratch/ambrad/test/test.kvf'], ncpu);
  pr('%s\n', estr);
  doex = N <= 3000;
  c = WriteKvf(N, 1e-5, 2); tic(); system(estr); toc
  if (doex) Ba = Extract(c.write_hmat_filename); fi(1); iml(Ba); cb; end
  fi(2); ShowBc(c);
  %return
  id = hm_mvp('i', c.write_hmat_filename);
  pr('true %e\n', sqrt(hm_mvp('fronorm2', id)));
  hm_mvp('c', id);
  %return;
  otol = 1e-3;
  c = WriteKvf(N,    otol, 2, c); tic; system(estr); toc
  if (doex) B1 = Extract(c.write_hmat_filename); end
  c = WriteKvf(N, 10*otol, 2, c); tic; system(estr); toc
  if (doex) B2 = Extract(c.write_hmat_filename); end
  if (doex) pr('%e (want <= %e)\n', relerr(B1, B2), 10*otol); end
end

function [me ma] = test_aca()
  m = 300;
  n = 370;
  A = randn(m, n);
  [U S V] = svd(A, 'econ');
  s = diag(S);
  s = logspace(0, -15, length(s));
  A = U*diag(s)*V';
  
  function B = Bfn(rs, cs)
    B = A(rs,cs);
  end
  
  blk = [10 11 m-20 n-21];
  scale = -1;
  use_rel_tol = 0;
  err = 1e-3;
  t = tic();
  [me.U me.V] = mextaca('aca', blk, @Bfn, scale, use_rel_tol, err);
  etme = toc(t);
  t = tic();
  [ma.U ma.V] = hm('BebACA', @(rs,cs) Bfn(rs + blk(1) - 1, cs + blk(2) - 1),...
                   scale, blk(3), blk(4), ~use_rel_tol, err);
  etma = toc(t);
  B = Bfn(blk(1):blk(1)+blk(3)-1, blk(2):blk(2)+blk(4)-1);
  test_cmp(ma, me, B, etma, etme);
  return;
  
  err = err/10;
  t = tic();
  [me.U me.V] = mextaca('aca', blk, @Bfn, scale, use_rel_tol, err, me.U, me.V);
  etme = toc(t);
  t = tic();
  [ma.U ma.V] = hm('WarmStartBebACA',...
                   @(rs,cs) Bfn(rs + blk(1) - 1, cs + blk(2) - 1),...
                   blk(3), blk(4), ~use_rel_tol, err, ma.U, ma.V);
  etma = toc(t);
  B = Bfn(blk(1):blk(1)+blk(3)-1, blk(2):blk(2)+blk(4)-1);
  dU = nan; dV = nan;
  if (size(ma.U,2) == size(me.U,2)) dU = relerr(ma.U, me.U); end
  if (size(ma.V,2) == size(me.V,2)) dV = relerr(ma.V, me.V); end
  fprintf(1, 'B vs U V'': %e %e | U, V vs U, V: %e %e\n',...
          relerr(B, ma.U*ma.V'), relerr(B, me.U*me.V'),...
          dU, dV);
  fprintf(1, 'tot et mat/mex: %1.2f\n', etma/etme);
  fprintf(1, 'rank: %d of %d\n', size(ma.U,2), min(m,n));

  err = 100*err;
  t = tic();
  [me.U me.V] = mextaca('compressqr', me.U, me.V, use_rel_tol, err);
  etme = etme + toc(t);
  t = tic();
  [ma.U ma.V] = CompressQR(ma.U, ma.V, err, ~use_rel_tol);
  etma = etma + toc(t);
  test_cmp(ma, me, B, etma, etme);
end

function [me ma] = test_lra()
  m = 300;
  n = 371;
  A = randn(m, n);
  [U S V] = svd(A, 'econ');
  s = diag(S);
  s = logspace(0, -15, length(s));
  A = U*diag(s)*V';
  
  function B = Bfn(rs, cs)
    B = A(rs,cs);
  end
  
  blk = [10 11 m-20 n-21];
  cqr = true;
  
  realp = 2;
  tol = 1e-3;
  t = tic();
  [me.B me.U me.V] = mextaca('lra', blk, @Bfn, tol, realp);
  etme = toc(t);
  rs = blk(1):blk(1)+blk(3)-1;
  cs = blk(2):blk(2)+blk(4)-1;
  opts.realp = realp;
  opts.bebaca_factor = 0.1;
  t = tic();
  [ma.U ma.V ma.B] = hm('LowRankApprox', @Bfn, rs, cs, tol, 'bebaca', true,...
                        cqr, opts);
  etma = toc(t);
  B = Bfn(rs, cs);
  test_cmp(ma, me, B, etma, etme);
  
  realp = 2;
  otol = tol;
  oma = ma;
  tol = tol*10;
  t = tic();
  [me.B me.U me.V] = mextaca('lra', blk, @Bfn, tol, realp, me.U, me.V, otol);
  etme = etme + toc(t);
  t = tic();
  opts.oh_in.tol = otol;
  opts.new_tol = tol;
  opts.ob = struct('B',[],'U',ma.U,'V',ma.V);
  [ma.U ma.V ma.B] = hm('LowRankApprox', @Bfn, rs, cs, tol, 'bebaca', true,...
                        cqr, opts);
  etma = etma + toc(t);
  test_cmp(ma, me, B, etma, etme);
  
  %return;
  [me.B me.U me.V] = mextaca('lra', blk, @Bfn, tol, realp, oma.U, oma.V, otol);
  test_cmp(ma, me, B, etma, etme);
end

function test_cmp(ma, me, B, etma, etme)
  dU = nan; dV = nan;
  if (size(ma.U,2) == size(me.U,2)) dU = relerr(abs(ma.U), abs(me.U)); end
  if (size(ma.V,2) == size(me.V,2)) dV = relerr(abs(ma.V), abs(me.V)); end
  daB = nan; deB = nan;
  if (~isempty(ma.U)) daB = relerr(B, ma.U*ma.V'); end
  if (~isempty(me.U)) deB = relerr(B, me.U*me.V'); end
  fprintf(1, 'B vs U V'': %e %e | U, V vs U, V: %e %e\n', daB, deB, dU, dV);
  fprintf(1, 'tot et mat/mex: %1.2f\n', etma/etme);
  fprintf(1, 'rank: %d %d of %d\n', size(ma.U,2), size(me.U,2), min(size(B)));
end

function test_BfroFromNearDiag(fn)
  id = hm_mvp('i', fn);
  nf2t = hm_mvp('fronorm2', id);
  hm_mvp('c', id);
  fprintf(1,'-> %e\n', sqrt(nf2t));
  
  fid = fopen(fn, 'r');
  [p q m n nb realp] = hm('ReadHeader', fid);
  nf2 = 0;
  pnf2 = 0;
  max_sz = 0;
  for (i = 1:nb)
    [b nnz] = hm('ReadBlock', fid, realp);
    rs = b.r0 + (0:b.m-1);
    cs = b.c0 + (0:b.n-1);
    if (~isempty(intersect(rs, cs)))
      if (b.m*b.n > max_sz) max_sz = b.m*b.n; end
      nf2 = nf2 + sum(b.B(:).^2);
      if (pnf2/nf2 < 0.1)
        pnf2 = nf2;
        fprintf(1,'%1.1e (%d) ', sqrt(nf2), max_sz);
      end
    end
  end
  fclose(fid);
  fprintf(1, 'max_sz = %d\n', max_sz);

  fprintf(1, '%e %e (%1.1f)\n', sqrt(nf2t), sqrt(nf2), 100*sqrt(nf2/nf2t));
end

function j = SvdChooseRank(s, rerr, rm_mrem)
  if (~rm_mrem)
    j = find(sqrt(cumsum(s(end:-1:1).^2))/norm(s,2) <= rerr, 1,'last');
  else
    j = find(sqrt(cumsum(s(end:-1:1).^2)) <= rerr, 1,'last');
  end
  if (isempty(j)) j = 0; end
  j = length(s) - j;
  % In practice I prefer not to have rank-0 blocks.
  if (j == 0) j = 1; end
end
  
function [U V] = CompressQR(U,V,rerr,rm_mrem)
  n = size(U,2);
  if (n > 1000) return; end
  oU = U; oV = V;
  [U UR] = qr(U,0);
  [V VR] = qr(V,0);
  [U1 s V1] = svd(UR*VR');
  s = diag(s);
  j = SvdChooseRank(s,rerr,rm_mrem);
  if (j == length(s))
    U = oU;
    V = oV;
    return;
  end
  U1 = U1(:,1:j)*diag(s(1:j));
  V1 = V1(:,1:j);
  U = U*U1;
  V = V*V1;
end

function nf2 = NormFro2UV(U, V)
  if (false)
    nf2 = sum(sum(V'.*((U'*U)*V')));
  else
    nf2 = 0;
    for (i = 1:size(U,2))
      nf2 = UpdateNormFro2UV(nf2, U(:,1:i-1), V(:,1:i-1), U(:,i), V(:,i));
    end
  end
end

function nf2 = UpdateNormFro2UV(nf2, U, V, u, v)
% u and v are the new columns that will be appended to U and V.
  nf2 = nf2 + 2*(u'*U)*(V'*v) + dot(u,u)*dot(v,v);
end

function CmpAcc(fns)
% Just for fun, compare H-matrices of varoius accuracies with a matrix in
% which the permuted form (with energy squeezed toward the diag) is reduced
% to a band such that the total storage is the same as the lowest-res
% H-matrix.
  
  % Read in the H-matrices.
  for (i = 1:numel(fns))
    in = hm('HmatInfo', fns{i});
    if (in.m*in.n > 4000^2) error('too big'); end
    [id nnz] = hm_mvp('i', fns{i});
    Bs{i} = hm_mvp('extract', id, 1:in.m, 1:in.n);
    cfs(i) = in.m*in.n/nnz;
    hm_mvp('c', id);
    pr('cf: %1.1f\n', cfs(i));
  end

  % Get the banded matrix. Physically, this is like completely ignoring
  % long-distance interactions.
  B = Bs{1};
  Bi = B(in.p, in.q);
  Bi = MakeBanded(Bi, round(size(B,1)/(2*max(cfs))));
  Bb(in.p, in.q) = Bi;
  Bs{end+1} = Bb;
  clear B Bi Bb;

  % Image of the banded matrix permuted back to physical space.
  fi(1); iml(Bs{end});
  % Select the middle column and compare the profiles.
  k = round(size(Bs{1}, 1)/2);
  fi(2); clf;
  for (i = 1:numel(Bs))
    b = Bs{i}(:,k);
    % Print 2-norm, 1-norm, and sum. The 2-norms tend to be about the same,
    % but the 1-norm and sum of the banded matrix is generally obviously
    % smaller than those of the H-matrices.
    pr('%e %e %e\n', norm(b, 2), norm(b, 1), sum(b));
    plot(b); hold all;
  end
end

function B = MakeBanded(B, k)
  for (i = 1:size(B,2))
    B(1:i-k-1,i) = 0;
    B(i+k+1:end,i) = 0;
  end
end

function [c d] = test_Elastostatics()
  c.x = sort(randn(1,10));
  c.eta = sort(randn(1,5));
  c.depth_min = 0;
  c.y_min = 0;
  c.dipdeg = 12;
  kvf('Write', '/scratch/ambrad/test/esin.kvf', c, 1);
  system(['./test_Elastostatics /scratch/ambrad/test/esin.kvf ',...
          '/scratch/ambrad/test/esout.kvf']);
  d = kvf('Read', '/scratch/ambrad/test/esout.kvf');
  
  c.eta = c.eta - c.eta(1);
  [X Eta] = meshgrid(c.x, c.eta);
  Y = cosd(c.dipdeg)*Eta;
  Z = sind(c.dipdeg)*Eta; Z = Z - max(Z(:));
  clf; hold all;
  plot3(X, Y, Z, 'r.');
  plot3(d.ctrs(1,:), d.ctrs(2,:), d.ctrs(3,:), 'b.');
  h = quiver3(d.ctrs(1,:), d.ctrs(2,:), d.ctrs(3,:),...
              d.adips(1,:), d.adips(2,:), d.adips(3,:));
  set(h, 'color', 'b');
  h = quiver3(d.ctrs(1,:), d.ctrs(2,:), d.ctrs(3,:),...
              d.astrikes(1,:), d.astrikes(2,:), d.astrikes(3,:));
  set(h, 'color', 'r');
  h = quiver3(d.ctrs(1,:), d.ctrs(2,:), d.ctrs(3,:),...
              d.normals(1,:), d.normals(2,:), d.normals(3,:));
  set(h, 'color', 'g');
  xlabel('x'); ylabel('y'); zlabel('z');
  hold off;
  axis equal;
end

function gf = StudyNrepeats()
  system('./test_nrepeats');
  gf = MatIo('ReadMat', 'nr_gf.mat');
  rs = [2 4];
  sp(211);
  nr = (size(gf,2) - 1)/2;
  ctr = nr + 1;
  plot(gf(1,:), log10(abs(mop(gf(rs,:),'/',gf(rs,ctr)))), '.-',...
       gf(1,:), log10(abs(1./(gf(1,:)/gf(1,ctr)).^3)), 'ko--');
  legend({'along-dip' 'normal' '1/r^3'}, 'location', 'ne');
  at; grid on;
  sp(212);
  cs = zeros(numel(rs), nr+1);
  for (i = 0:nr)
    cs(:,i+1) = sum(gf(rs, ctr-i:ctr+i), 2);
  end
  plot(0:nr, log10(abs(mop(mop(cs,'-',cs(:,end)), '/',...
                           cs(:,end)))), '.-');
  at; grid on;
end

function StudyNrepeatsH(N, tol, ncpu, nrepeat)
  if (N > 3000) error('N is to large'); end
  prob = 3;
  estr = sprintf(['mpirun -np %d ./bin/hmmvpCompress ',...
                  '/scratch/ambrad/test/test.kvf'], ncpu);
  pr('%s\n', estr);

  c = WriteKvf(N, tol, prob, [], 'nrepeat', 10); tic(); system(estr); o.et1 = toc();
  pr('et: %f\n', o.et1);
  o.B1 = Extract(c.write_hmat_filename);
  o.cf1 = GetCompressFactor(c.write_hmat_filename);
  c = WriteKvf(N, tol, prob, [], 'nrepeat', nrepeat); tic(); system(estr); o.et2 = toc();
  pr('et: %f\n', o.et2);
  o.B2 = Extract(c.write_hmat_filename);
  o.cf2 = GetCompressFactor(c.write_hmat_filename);
  
  o.c = c;
  sp(221); iml(o.B1); cb; c = caxis();
  sp(222); iml(o.B2); cb; caxis(c);
  sp(223); iml((o.B1-o.B2)./o.B1); cb;
  pr('nwre = %e\n', relerr(o.B1, o.B2));

  pr('%f vs %f for speed factor %1.1f\n', o.et1, o.et2, o.et1/o.et2);
end

function RunNrH(N, tol, ncpu, nrepeat)
  prob = 3;
  estr = sprintf(['mpirun -np %d ./bin/hmmvpCompress ',...
                  '/scratch/ambrad/test/test.kvf'], ncpu);
  pr('%s\n', estr);
  c = WriteKvf(N, tol, prob, [], 'nrepeat', nrepeat);
  tic(); system(estr); o.et1 = toc();
  pr('et: %f\n', o.et1);
end

function A = StudyBlocks(fn)
  [bs p q nnz m n] = hm('ReadBlocks', fn);
  A = zeros(m, n);
  for (i = 1:len(bs))
    A(bs(i).r0 + (1:bs(i).m), bs(i).c0 + (1:bs(i).n)) = rand;
  end
  ip(p+1) = 1:m;
  iq(q+1) = 1:n;
  imagesc(A(ip,iq));
end

function A = StudyDiag(fn)
  function id = IsDiagBlock(b, p, q)
    if (false)
      r = p(b.r0 + (1:b.m));
      c = q(b.c0 + (1:b.n));
      id = numel(r) == numel(c) && all(r == c);
    else
      sz = 3;
      id = b.m <= sz && b.n <= sz;
    end
  end
  
  [bs p q nnz m n] = hm('ReadBlocks', fn);
  A = zeros(m, n);
  for (i = 1:len(bs))
    if (IsDiagBlock(bs(i), p, q))
      A(bs(i).r0 + (1:bs(i).m), bs(i).c0 + (1:bs(i).n)) = 1;
    end
  end
  ip(p+1) = 1:m;
  iq(q+1) = 1:n;
  imagesc(A);
end

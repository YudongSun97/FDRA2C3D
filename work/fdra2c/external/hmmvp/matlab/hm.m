function varargout = hm(varargin)
% hm('BlackBoxCompress',bs,permp,permq,Bfn,rerr,savefn)
%   Compress the blocks in the H-matrix approximating B. All options are
% chosen to be as good as I think they can be, and 'rerr' has a very obvious
% meaning. Call this routine rather than hm('Compress') unless you want to
% play around with options.
%
% INPUTS
%   bs,p,q are the set of blocks and the permutation obtained from
%     [bs p q] = hmcc_hd(D,R);
%   Bfn is a function having signature
%     B = Bfn(rs,cs),
% where B is the BEM matrix; Bfn can also simply be the matrix when it is not
% too large. rs and cs are row and column indices.
%   rerr is the requested relative error. It's the bound on the maximum
% relative normwise relative error in the matrix-vector product. It should be
% set to something like 1e-2 to 1e-4, with 1e-3 being about right most of the
% time.
%   savefn is the file to which to save the H-matrix data.
%
%   In more detail, the H-matrix is created so that the error between the
% true matrix and its approximation is such that the error in the
% matrix-vector product is bounded as
%     ||E x||_2 <= rerr ||B x||_2.
% This is just the normwise relative error of the matrix product.
%
% ------------------------------------------------------------------------------
%
% hm('Compress',bs,p,q,Bfn,rerr,savefn)
%   Compress the blocks in the H-matrix approximating B.
%
% INPUTS
%   The required arguments are the same as for BlackBoxCompress except:
%     rerr is the requested relative error. The meaning of rerr depends on
%     options; see below.
%
% OPTIONAL
%   Optional arguments are passed as a string-value pair.
%   'old_hmat_fn': Provide an H-matrix file from a previous call to
% hm('Compress') to speed up constructing the new H-matrix. bs and perms must
% be identical for the two H-matrices.
%   'tol_method': This option defines 'rerr'. Values are 'Bfro' [defualt],
% 'BinvFro', and 'abs'. If 'BinvFro' is used, either of the options 'BinvFro'
% or 'old_hmat_fn' must be provided.
%   'Bfro': Provide the (estimated) Frobenius norm of B. This option is used
% with 'tol_method' set to 'Bfro'.
%   'BinvFro': Provide the (estimated) Frobenius norm of inv(B). This method
% is used with 'tol_method' set to 'BinvFro'.
%
%   'rerr' has the following interpretation. If 'tol_method' is 'Bfro', then
% the error E between the true matrix B and its approximation is such that
% the error in the matrix-vector product is bounded as
%     ||E x||_2 <= rerr ||B||_F ||x||_2.
%   If 'tol_method' is 'BinvFro', then the error is bounded as
%     ||E x||_2 <= rerr ||B x||_2,
% which is a tighter bound than the first one on the relative error.
%   If 'tol_method' is 'abs', then the error is bounded as
%     ||E x||_2 <= rerr ||x||_2.
% This option is appropriate only if the maximum allowable absolute error in the
% MVP given a unit vector x is known.
%
% ----
%
% Bfro = hm('EstBfroLowBndCols',Bfn,p,q,rs,cs)
%   Estimate a lower bound on norm(Bfn,'fro') by using elements centered on the
% main diagonal in each column of the permuted matrix. The permuted matrix,
% where the permutation is given by hmcc_hd, tends to have most of its energy
% near the main diagonal.
%
% INPUTS
%   Bfn(r,c) gives the rows and cols r,c.
%   p,q are the permutation matrices returned by hmcc_hd.
%   rs,cs are the indices containing the diagonal elements in the unpermuted
%     B. For example, if B is structurally symmetric, then rs = cs = 1:m.
%
% OUTPUTS
%   Bfro: Estimated lower bound on ||B||_F.
%
% ----
%
% info = hm('HmatInfo',fn)
%   Return matrix size and row and col permutations of the H-matrix stored in
% the file fn, as well as other data.
%
% ----
% 
% fn = hm('FroNorm',fn)
%   Return the Frobenius norm of the H-matrix stored in the file fn.
%
% ----
%
%   The following methods perform some higher-level operations with
% H-matrices. In some problems, a client may construct a matrix using
% multiple H-matrices, for example, as a result of geometric
% symmetries. Consequently, the call hm_mvp('mvp',id,x) may not be
% sufficient; rather, the client needs to implement the matrix-vector
% product.
%   To allow this, a function handle can be passed in an 'opts' struct that
% has the same signature as 'hm_mvp'. Set opts.hm_mvp to the function handle
% to use this functionality.
%
% f = hm('BlocksPrecon',id,opts)
%   Build a preconditioner based on the full blocks. f has fields 'L' and
%   'U', which are the factors returned by ilu.
%
% [x o] = SolveBxb(id,b,reltol,opts)
%   Solve B x = b using GMRES. reltol is a relative error tolerance. opts is
% struct having the following optional fields:
%     .hm_mvp: As described above.
%     .is: Use B(is,is) rather than the full B.
%     .precon: Provide a preconditioner. It can have one of three forms:
%        1. A struct having fields L and U, e.g., from 'BlocksPrecon'.
%        2. A struct having fields L, U, P, Q, R, presumably from
%               [L U P Q R] = lu(...).
%        3. A function handle for a function having signature
%               x = @(b).
%     .dir_mvp: 'mvp' [default] to solve B x = b; 'tmvp' to solve B' x = b.
%     .x0: Starting point for GMRES [default is 0s].
%     .restart: 'restart' parameter for GMRES. Default is 20. Probably worth
%        changing only if the memory to hold 20 mesh-sized vectors is a
%        concern.
%
% ce = CondEst1(id,opts)
%   Like Matlab's condest1.
%
% [ce re] = CondEstFro(id,opts)
%   Estimate cond(B,'fro'). opts has the following optional fields:
%     .hm_mvp: As described above.
%     .reltol: the desired error level. It should be set no smaller than
%       about 0.01 [default].
%     .max_nits: the maximum allowed number of iterations.
%     .precon: 0 for no preconditioning; 1 [default] to request a
%       preconditioner; or a preconditioner struct to pass to SolveBxb.
%   ce is the estimate.
%   re is the jackknife estimate of the relative error in ce. It probably
% underestimates the error by a bit, but generally that's not a problem.
%
% [ce re] = NormInvEstFro(id,opts)
%   Estimate norm(inv(B),'fro'). Inputs and 're' are the same as for
% CondEstFro.

% Andrew M. Bradley ambrad@cs.stanford.edu
% CDFM, Geophysics, Stanford
  [varargout{1:nargout}] = feval(varargin{:});

% ------------------------------------------------------------------------------
% Block compression

function BlackBoxCompress(bs,permp,permq,Bfn,rerr,savefn)
  t = tic();
  n = length(permp);
  eBfro = hm('EstBfroLowBndCols',Bfn,permp,permq,1:n,1:n);
  savefn_tmp = [savefn '.tmp'];
  hm('Compress',bs,permp,permq,Bfn,max(1e-5*rerr,1e-8),savefn_tmp,'Bfro',eBfro);
  [id hnnz] = hm_mvp('init',savefn_tmp);
  in = hm('HmatInfo',savefn_tmp);
  cf = in.m*in.n/hnnz;
  Bfro = sqrt(hm_mvp('fronorm2',id));
  [ce re] = hm('CondEstFro',id,struct('reltol',0.02,'max_nits',50));
  hm_mvp('cleanup',id);
  et = toc(t);
  fprintf(1,['Evaluation H-matrix:\n',...
	     '  time: %1.1f\n',...
	     '  compression factor: %1.1f\n',...
	     '  ||B||_F: %1.1e (rough) %1.1e (refined)\n',...
	     '  cond(B,''fro''): %1.1e (+/- %1.1f%%)\n'],...
	  et,cf,eBfro,Bfro,ce,100*re);
  
  t = tic();
  hm('Compress',bs,permp,permq,Bfn,rerr,savefn,...
     'old_hmat_fn',savefn_tmp,'tol_method','BinvFro','BinvFro',ce/Bfro);
  [id hnnz] = hm_mvp('init',savefn);
  hm_mvp('cleanup',id);
  cf = in.m*in.n/hnnz;
  fprintf(1,['Final H-matrix:\n',...
	     '  time: %1.1f\n',...
	     '  compression factor: %1.1f\n'],...
	  toc(t),cf);

function tol_d_rerr = Compress(bs,permp,permq,Bfn,rerr,savefn,varargin)
%   precision: 1 (single), [2] (double)
  [method verbose rerr_method cqr bebaca_factor realp old_hmat_fn,...
   tol_method Bfro BinvFro] =...
      process_options(...
	  varargin,'method','bebaca','verbose',1,'rerr_method','mrem',...
	  'compress_qr',1,'bebaca_factor',0.1,'precision',1,'old_hmat_fn','',...
	  'tol_method','Bfro','Bfro',[],'BinvFro',[]);
  method = lower(method);
  rm_mrem = strcmp(rerr_method,'mrem') || strcmp(rerr_method,'global');
  have_old_hmat = ~isempty(old_hmat_fn);
  
  % Get dims of B for progress indicator.
  Bnr = length(permp); Bnc = length(permq);
  Bn = Bnr*Bnc;

  % Do what we can with the old H-matrix.
  if (have_old_hmat)
    if (strcmp(old_hmat_fn,savefn))
      error('New and old hmat filenames should not be the same.');
    end
    try
      oh_in = hm('HmatInfo',old_hmat_fn);
    catch
      error(sprintf('%s is not a valid H-matrix file.',old_hmat_fn));
    end
    if (oh_in.m ~= Bnr || oh_in.n ~= Bnc || oh_in.nb ~= length(bs) ||...
	~all(permp(:) == oh_in.p(:)) || ~all(permq(:) == oh_in.q(:)) ||...
	oh_in.rerr_method ~= rm_mrem)
      oh_err = sprintf(...
	  'H-matrix file %s is not compatible with this H-matrix.',...
	  old_hmat_fn);
      error(oh_err);
    end
    tol_d_rerr = 1;
    if (rm_mrem)
      switch (tol_method)
       case 'Bfro'
	if (isempty(Bfro))
	  % Get an accurate ||B||_F.
	  id = hm_mvp('init',old_hmat_fn);
	  Bfro = sqrt(hm_mvp('fronorm2',id));
	  hm_mvp('cleanup',id);
	end
	tol_d_rerr = Bfro;
       case 'BinvFro'
	if (isempty(BinvFro))
	  id = hm_mvp('init',old_hmat_fn);
	  BinvFro = NormInvEstFro(id);
	  hm_mvp('cleanup',id);
	end
	tol_d_rerr = 1/BinvFro;
       case 'abs'
        tol_d_rerr = 1;
       otherwise
        error(sprintf('%s is not a value for tol_method.',tol_method));
      end
    end
    if (rerr*tol_d_rerr == oh_in.tol)
      error(sprintf('Requested tolerance is same as for %s.\n',old_hmat_fn));
    end
  end
  
  % Figure out the user tolerance method.
  if (rm_mrem)
    switch (tol_method)
     case 'Bfro'
      if (isempty(Bfro))
	t = tic;
	% Estimating norm(B,'fro') is slow, but it's worth it.
	[Bfro sigma2] = EstBfroCols(...
	    @(cs)Bfn(1:Bnr,cs),Bnc,max(3000,ceil(0.03*sqrt(Bn))),1/50);
	fprintf(1,'et EstBfroCols = %f  Bfro = %1.18e  sqrt(sigma2) = %e\n',...
		toc(t),Bfro,sqrt(sigma2));
	% To prevent catastrophic errors, reset sigma2 if it's not small
	% enough. We likely won't achieve the user-requested tolerance, but at
	% least we won't bomb here. This should essentially never happen unless
	% B is really oddly behaved -- among other things, not at all smooth.
	if (2*sqrt(sigma2) > Bfro) sigma2 = 0; end
	% To be conservative, 2 std(estimator) from the estimator.
	Bfro = Bfro - 2*sqrt(sigma2);
      end
      tol_d_rerr = Bfro;
     case 'BinvFro'
      if (isempty(BinvFro))
	% We can't estimate ||inv(B)||_F without already having an H-matrix.
	error(['If rerr_method == mrem and tol_method == BinvFro, then\n',...
	       'BinvFro must be provided.\n']);
      end
      tol_d_rerr = 1/BinvFro;
     case 'abs'
      tol_d_rerr = 1;
     otherwise
      error(sprintf('%s is not a value for tol_method.',tol_method));
    end
    % We want tol/N, where B is NxN.
    tol_divN = rerr*tol_d_rerr/sqrt(Bn);
    if (tol_divN < 10*sqrt(eps)) realp = 2; end
  else
    tol_d_rerr = 1;
    % With BREM, it's not quite clear how to set this since we're not using what
    % is essentially an abs tol. I'll just go with this for now:
    if (rerr < 10*sqrt(eps)) realp = 2; end
  end
  
  savefn = Filename(savefn);
  fid = fopen(savefn,'wb');
  if (fid == -1)
    error(sprintf('Compress: Can''t open %s for writing.\n',savefn));
  end
  WriteHeader(fid,permp,permq,realp,length(bs),rm_mrem,rerr*tol_d_rerr);
  
  opts.bebaca_factor = bebaca_factor;
  opts.realp = realp;
  
  p = 0;
  fprintf(1,'Compress progress: ');
  if (have_old_hmat)
    opts.oh_in = oh_in;
    opts.new_tol = rerr*tol_d_rerr;
    clear oh_in;
    ofid = fopen(old_hmat_fn,'rb');
    ReadHeader(ofid);
    Cast = RealTypeFn(opts.realp);
  end
  for (i = 1:length(bs))
    % Progress
    p1 = floor(100*i/length(bs));
    if (verbose >= 1 && p1 > p)
      fprintf(1,'%d ',p1);
      p = p1;
    end

    if (false)
      nr = bs(i).m;
      nc = bs(i).n;
      r = permp(bs(i).r0 + (0:nr-1));
      c = permq(bs(i).c0 + (0:nc-1));
    else
      nr = bs(3,i);
      nc = bs(4,i);
      r = permp(bs(1,i) + (0:nr-1));
      c = permq(bs(2,i) + (0:nc-1));
    end
    
    if (have_old_hmat)
      opts.ob = ReadBlock(ofid,opts.oh_in.realp);
      % Cast to current type in case the old H-mat's type is different.
      flds = {'B' 'U' 'V'};
      for (j = 1:length(flds)) opts.ob.(flds{j}) = Cast(opts.ob.(flds{j})); end
    end
    
    if (nr == 1 || nc == 1)
      if (have_old_hmat) B = opts.ob.B;
      else               B = Bfn(r,c);  end
      U = []; V = [];
    else
      if (~rm_mrem) rerr1 = rerr;
      else          rerr1 = sqrt(nr*nc)*tol_divN; end
      [U V B err] = LowRankApprox(Bfn,r,c,rerr1,method,rm_mrem,cqr,opts);
    end

    WriteBlock(fid,realp,bs(1,i),nr,bs(2,i),nc,U,V',B);
        
    if (verbose >= 2 && nr*nc > 10000)
      fprintf(1,'%dx%d %3.1f %3.1e\n',nr,nc,...
	      nr*nc/(prod(size(U)) + prod(size(V)) + prod(size(B))),err);
    end
  end
  fprintf(1,'\n');
  if (have_old_hmat) fclose(ofid); end
  fclose(fid);
  
function [U V B err] = LowRankApprox(Bfn,r,c,rerr,method,rm_mrem,cqr,opts)
% If rm_mrem, then rerr is an absolute error tol.
  if (nargin < 8) opts = []; end
  if (~isfield(opts,'realp')) opts.realp = 2; end
  Cast = RealTypeFn(opts.realp);
  nr = length(r);
  nc = length(c);
  err = inf;
  B = [];

  hoh = isfield(opts,'ob');
  houv = false;
  if (hoh)
    old_tol = opts.oh_in.tol;
    new_tol = opts.new_tol;
    % In three of four cases, we don't need to evaluate Bfn. In two of those,
    % we can do a quick calculation and return:
    if (new_tol < old_tol && ~isempty(opts.ob.B))
      B = Cast(opts.ob.B);
      err = 0; U = []; V = [];
      return;
    elseif (new_tol > old_tol && isempty(opts.ob.B))
      if (size(opts.ob.U,2) > 1)
	alpha = old_tol/new_tol;
	beta = 1 - alpha;
	if (~rm_mrem) beta = beta / (1 + alpha*rerr); end	
	[U V] = CompressQR(opts.ob.U,opts.ob.V,beta*rerr,rm_mrem);
      else
	U = opts.ob.U; V = opts.ob.V;
      end
      B = [];
      return;
    end
    % If we're still here, either
    %   1. we need greater accuracy and have only ob.U, ob.V or
    %   2. we need less accuracy and have ob.B.
    % Case 1 requires modifying the ACA routine. Flag that here:
    houv = new_tol < old_tol && isempty(opts.ob.B);
    % Case 2 is easy: in what follows, simply provide B instead of Bfn:
    if (~houv) B = opts.ob.B; end
  end
  
  switch (method)
   
   case 'svd'
    if (isempty(B)) B = Cast(Bfn(r,c)); end
    [U S V] = svd(B,'econ');
    s = diag(S);
    j = SvdChooseRank(s,rerr,rm_mrem);
    if ((nr + nc)*j < nr*nc)
      U = U(:,1:j)*diag(s(1:j));
      V = V(:,1:j);
      err = norm(B - U*V','fro');
      if (~rm_mrem) err = err/norm(B,'fro'); end
      B = [];
    else
      U = []; V = [];
    end
       
   case 'bebaca'
    if (cqr) alpha = 0.5; else alpha = 1; end
    if (max(nr,nc) <= 250)
      [U V B err] = LowRankApprox(Bfn,r,c,rerr,'svd',rm_mrem,cqr,opts);
    else
      if (isempty(B)) Bfn1 = @(i,j) Cast(Bfn(r(i),c(j)));
      else            Bfn1 = Cast(B); end
      % We set this first error to ~1/10 of the requested error because ACA
      % sometimes is not accurate enough by ~<=10x. CompressQR mostly makes
      % up for this. (bebaca_factor is nominally 1/10, but can have another
      % value.)
      rerr1 = rerr;
      if (cqr) rerr1 = rerr*opts.bebaca_factor; end
      if (~houv)
        %todo Bscale is not important, right?
	Bscale = -1;
	[U V flag] = BebACA(Bfn1,Cast(Bscale),nr,nc,rm_mrem,alpha*rerr1);
	% We used to fail on ~flag. But the correct interpretation of ~flag
        % is that the matrix is exactly low rank and we found that rank. So
        % ignore it unless we think of a reason not to.
	flag = true;
      else
	[U V flag] = WarmStartBebACA(...
	    Bfn1,nr,nc,rm_mrem,alpha*rerr1,opts.ob.U,opts.ob.V);
      end
      if (~flag || prod(size(U)) + prod(size(V)) > nr*nc)
	B = Bfn(r,c);
	err = 0; U = []; V = [];
      else
	B = [];
	if (cqr && ~isempty(U))
	  beta = 1 - alpha;
	  if (~rm_mrem) beta = beta / (1 + alpha*rerr); end
	  [U V] = CompressQR(U,V,beta*rerr,rm_mrem);
	end
      end
    end

  end
  
function j = SvdChooseRank(s,rerr,rm_mrem)
  if (~rm_mrem)
    j = find(sqrt(cumsum(s(end:-1:1).^2))/norm(s,2) <= rerr, 1,'last');
  else
    j = find(sqrt(cumsum(s(end:-1:1).^2)) <= rerr, 1,'last');
  end
  if (isempty(j)) j = 0; end
  j = length(s) - j;
  % In practice I prefer not to have rank-0 blocks.
  if (j == 0) j = 1; end
  
function [U V] = CompressQR(U,V,rerr,rm_mrem)
  global use_mex;
  if (use_mex)
    [U V] = mextaca('compressqr', double(U), double(V), ~rm_mrem,...
                    double(rerr));
    return;
  end
  
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
  
function [U V flag] = BebACA(Bfn,scale,M,N,rm_mrem,err,U,V)
% Essentially a duplicate of ACAr in Bebendorf's AHMED/ACA.h.
%   To get single-precision output, Bfn and err should both be 'single'.

% TODO: I could write a separate file storing Z arrays at termination. Then
% the warm-started version wouldn't waste time trying rows already
% used. Experiments show that isn't a huge issue, though.

  global use_mex;
  if (use_mex)
    flag = 1;
    blk = [1 1 M N];
    Bfn1 = @(i,j) double(Bfn(i,j));
    if (nargin < 7)
      U = [];
      V = [];
    end
    [U V] = mextaca(...
        'aca', blk, Bfn1, double(scale), ~rm_mrem, double(err),...
        double(U), double(V));
    return;
  end

  flag = 0;
  if (nargin < 7)
    U = ones(0,0,class(err)); % Inherit precision of 'err'
    V = U;
    k1 = 1;
  else
    k1 = size(U,2) + 1;
  end
  Z = zeros(M,1);
  if (isempty(U))
    irow = floor(M/2) + 1;
  else
    [~,irow] = max(abs(U(:,end)));
  end
  if (isa(err,'single')) scale_fac = 1e1*eps(single(1));
  else scale_fac = 1e2*eps(double(1)); end
  for (k = k1:min(M,N))
    while (true)
      % New row of R
      v = Bfn(irow,1:N);
      if (k > 1) v = v - U(irow,:)*V.'; end
      % Get column
      [~,icol] = max(abs(v));
      Z(irow) = nan;
      % Acceptable?
      if (abs(v(icol)) > scale_fac*scale) break; end
      % No, so try a new irow
      fnd = false;
      for (i = 1:M)
	irow = irow + 1; if (irow > M) irow = 1; end
	if (~isnan(Z(irow))) fnd = true; break; end
      end
      % Dealing with an exactly low-rank matrix.
      if (~fnd) return; end
    end
    % R(:,icol)
    u = Bfn(1:M,icol);
    if (k > 1) u = u - U*V(icol,:).'; end
    % Append u,v
    v = v/v(icol);
    U(:,end+1) = u;
    V(:,end+1) = v;
    % Next irow
    u1 = u;
    u1(isnan(Z)) = nan;
    [~,irow] = max(abs(u1));
    % Assess error
    Bfro2 = NormFro2UV(U,V);
    Rfro2 = NormFro2UV(u,v);
    % Good?
    if (rm_mrem)
      if (Rfro2 <= err^2) break; end
    else
      if (Rfro2 <= err^2*Bfro2) break; end
    end
    scale = sqrt(Rfro2/(M*N));
  end
  flag = 1;

function [U V flag] = WarmStartBebACA(Bfn,M,N,rm_mrem,err,U,V)
% Wrapper to BebACA when starting with U,V.
  scale = sqrt(NormFro2UV(U,V)/(M*N));
  [U V flag] = BebACA(Bfn,scale,M,N,rm_mrem,err,U,V);
  
function nf2 = NormFro2UV(U,V)
% Frobenius norm squared of U*V'.
  nf2 = sum(sum(V'.*((U'*U)*V')));
  
% ------------------------------------------------------------------------------
% Methods to estimate norm(B,'fro'). The ideal method to get a BREM H-matrix is
% as follows: Get a cheap initial estimate of Bfro. Compute a low-tol H-matrix
% using the estimate as the 'Bfro' optional argument. Use 'FroNorm' to get
% the norm of that H-matrix, and use it as a good estimate to compute the
% final H-matrix. If 'Bfro' is not passed to 'Compress', then 'EsBfroCols' is
% called with rather conservative, and so expensive, parameters.

function [Bfro sigma2 ncols] = EstBfroCols(Bfn,N,ncols0,max_ratio)
% Estimate norm(Bfn,'fro') by using a subset of the columns of Bfn.
% Input:
%   Bfn(c) is column c of B.
%   B has N columns.
%   ncols0: Start by using this many columns.
%   max_ratio: But don't stop until sqrt(sigma2)/Bfro <= max_ratio.
% Output:
%   Bfro: Estimate of ||B||_F.
%   sigma2: Variance of estimator.
%   ncols: Number of columns used.
  
  if (nargin < 5) max_ratio = 1; end
  
  pncols = 0;
  ncols = ncols0;
  idxs = [];
  Bcol2s = [];
  while (true)
    fprintf(1,'norm(B,''fro'') progress: ');
    plen = length(idxs);
    % Select a subset of columns randomly
    if (N < ncols)
      idxs = 1:N;
    else
      rp = randperm(N);
      rp = setdiff(rp,idxs);
      idxs = [idxs rp(1:(ncols - pncols))];
    end
    Ne = length(idxs) - plen;
    % Get the 2-norm (squared) of each column
    Bcol2s = [Bcol2s zeros(1,Ne)];
    p = 0;
    for (j = 1:Ne)
      p1 = floor(100*j/Ne); if (p1 > p) fprintf(1,'%d ',p1); p = p1; end
      i = idxs(plen + j);
      Bcol2s(plen + j) = sum(Bfn(i).^2);
    end
    fprintf(1,'\n');
    % Get mean and jackknife variance of this estimator
    EstFn = @(X) sqrt(N/length(X)*sum(X));
    [Bfro sigma2] = Jackknife(EstFn,Bcol2s);
    
    if (length(idxs) == N || sqrt(sigma2)/Bfro <= max_ratio) break; end
    pncols = ncols;
    ncols = 2*ncols;
  end

function Bfro = EstBfroLowBndCols(Bfn,p,q,rs,cs)
% Estimate a lower bound on norm(Bfn,'fro') by using elements centered on the
% main diagonal in each column of the permuted matrix. The permuted matrix,
% where the permutation is given by hmcc_hd, tends to have most of its energy
% near the main diagonal.
% Input:
%   Bfn(r,c) gives the rows and cols r,c.
%   p,q are the permutation matrices returned by hmcc_hd.
%   rs,cs are the indices containing the diagonal elements in the unpermuted
%     B. For example, if B is structurally symmetric, then rs = cs = 1:m.
% Output:
%   Bfro: Estimated lower bound on ||B||_F.
  
  % Rows of the main diagonal in the permuted matrix
  ip = invperm(p);
  rs = ip(rs);
  nr = length(p);
  nc = length(cs);
  
  % Figure out an appropriate n, the bandwidth
  N = 10;
  ns = zeros(1,N);
  prop = 0.99;
  Ks = unique(round(linspace(1,nc,N)));
  for (i = 1:N)
    k = Ks(i);
    b = Bfn(p,cs(k));
    ns(i) = GetAllEnergyAbove(b,rs(k),prop);
  end
  n = round(median(ns));
  %fprintf(1,'n = %d\n',n);
  
  Bfn1 = @(k) Bfn(p(max(1,rs(k)-n):min(nr,rs(k)+n)),cs(k));
  [Bfro sigma2 ncols] = EstBfroCols(Bfn1,nc,3000,1/50);
    
function [rs cs] = GetMainDiagIdxs(p,q,rs,cs)
% On entry, (rs,cs) are the rows and columns of elements of interest (eg, the
% main diagonal) in the original (unpermuted) matrix B. p,q are permutation
% matrices associated with the H-matrix and obtained from hmcc_hd. On exit,
% (rs,cs) are the rows and columns of the main diagonal in the permuted matrix.
  ip = invperm(p);
  iq = invperm(q);
  rs = ip(rs);
  cs = iq(cs);
  
function n = GetAllEnergyAbove(b,r,p)
% Find the minimum n such that for the indices
%     is = max(1,r-n):min(length(b),r+n),
% which are a window centered at index r,
%     sum(abs(b(is))) >= p*sum(abs(b)).
  % Finagle things in an efficient way so that...
  bcsf = cumsum(abs(b(r:end)));
  bcsb = [0; cumsum(abs(b(r-1:-1:1)))];
  bs = bcsf(end) + bcsb(end); % sum(abs(b))
  % ...cs is the cumulative sum defined by increasing n by 1 going outward
  % from r.
  cs = sum([bcsf(:)' bcsf(end)*ones(1,length(bcsb) - length(bcsf))
	    bcsb(:)' bcsb(end)*ones(1,length(bcsf) - length(bcsb))]);
  % Could do binary search on cs at this point, but we've already had to do
  % O(length(b)) work, so might as well use the built-in find.
  n = find(cs >= p*bs,1) - 1;
  
function [ev sigma2] = Jackknife(EstFn,X)
% X is KxN and contains N K-vector samples. EstFn is the estimator and has
% signature
%     y = EstFn(X),
% where X is Kx(N-1) and y is the estimate corresponding to this subsample.
%   ev is the estimate and sigma2 is the delete-1 jackknife variance of the
% estimator.
  Ne = size(X,2);
  if (false) % would use this if I add a 'vectorized' option
    % ith column of p is setdiff(1:Ne,i)
    i = 1:Ne;
    p = repmat(i',1,Ne);
    p(Ne*(i-1)+i) = [];
    p = reshape(p,Ne-1,Ne);
  end
  % Get estimate
  ev = EstFn(X);
  % Get estimates corresponding to each subsample
  Y = zeros(1,Ne);
  p = 1:Ne;
  for (i = 1:Ne)
    p1 = p;
    p1(i) = [];
    Y(i) = EstFn(X(:,p1));
  end
  % Jackknife variance of the estimator
  sigma2 = sum((Y - mean(Y)).^2)/Ne;
  
% ------------------------------------------------------------------------------
% Hmat operations

function [bs p sz] = SortBlocks(bs)
  sz = zeros(1,length(bs));
  for (i = 1:length(bs))
    sz(i) = bs(i).m*bs(i).n;
  end
  [~,p] = sort(-sz);
  sz = sz(p);
  bs = bs(p);

function y = Mvp(bs,p,q,x)
  p = p + 1;
  q = q + 1;
  x = x(q);
  y = zeros(size(p));
  for (i = 1:length(bs))
    b = bs(i);
    if (isempty(b.B))
      y(b.r0+(1:b.m)) = y(b.r0+(1:b.m)) + b.U*(x(b.c0+(1:b.n))'*b.V)';
    else
      y(b.r0+(1:b.m)) = y(b.r0+(1:b.m)) + b.B*x(b.c0+(1:b.n));
    end
  end
  y(p) = y;
  
function ip = invperm(p)
  ip = p;
  ip(p) = 1:length(p);
    
% ------------------------------------------------------------------------------
% Block I/O

function fn = Filename(fn)
  return; % decimal points in the filename mess this up
  [p n e] = fileparts(fn);
  if (isempty(e))
    fn = [fn '.dat'];
  end

function WriteHeader(fid,permp,permq,realp,nb,rm_mrem,tol)
  int = FilePrecision();
  write(fid, length(permp),   int);
  write(fid, length(permq),   int);
  write(fid, realp,           int);
  write(fid, permp - 1,       int);
  write(fid, permq - 1,       int);
  write(fid, nb,              int);
  write(fid, double(rm_mrem), int);
  % For MREM, In v0.8 to v0.13, we wrote double([rerr Bfro]), where Bfro is an
  % estimate. But now we write just tol, where ||dB|| <= tol. tol can be derived
  % in any way desired; it doesn't matter for purposes of recording data and
  % later recovering the tol. To indicate that this is >=v0.14, place 14.0 in
  % the second (now unused) double field. Perhaps we'll use it in a future
  % version.
  %  If BREM is used, then tol is the relative tolerance.
  write(fid, tol,             'double');
  write(fid, 14.0,            'double');
  
function WriteBlock(fid,realp,r0,m,c0,n,U,Vt,B)
  [int real char] = FilePrecision(realp);
  write(fid,   r0-1,      int );
  write(fid,   m,         int );
  write(fid,   c0-1,      int );
  write(fid,   n,         int );
  if (isempty(B))
    write(fid, 'U',       char);
    write(fid, size(U,2), int );
    write(fid, U(:),      real);
    write(fid, Vt(:),     real);
  else
    write(fid, 'B',       char);
    write(fid, B(:),      real);
  end
  
function [bs p q nnz m n] = ReadBlocks(savefn)
  fid = fopen(Filename(savefn),'rb');
  [p q m n nb realp] = ReadHeader(fid);
  [int real] = FilePrecision(realp);
  nnz = 0;
  for (i = 1:nb)
    [b bnnz] = ReadBlock(fid,realp);
    if (i == 1)
      % Allocate storage; otherwise, adding to bs slows this down by a lot.
      bt = b;
      bt.B = []; bt.U = []; bt.V = [];
      bs = repmat(bt,1,nb);
    end
    bs(i) = b;
    nnz = nnz + bnnz;
  end
  fclose(fid);
  
function [p q m n nb realp] = ReadHeader(fid)
  int = FilePrecision(); % FilePrecision establishes variable storage classes.
  m  = read(fid,1,int);  % Number of rows in the H-matrix.
  n  = read(fid,1,int);  % Number of cols.
  realp = read(fid,1,int); % Data are stored as float (1) or double (2).
  p  = read(fid,m,int);  % Row permutations from user to H-matrix indexing.
  q  = read(fid,n,int);  % Col permutations.
  nb = read(fid,1,int);  % Number of blocks.
  % Don't need these data: just eat them to move the file pointer forward.
  read(fid,1,int); read(fid,1,'double'); read(fid,1,'double');
  
function [b nnz] = ReadBlock(fid,realp)
  [int real char] = FilePrecision(realp);
  A = read(fid, 4, int);
  % The block covers rows b.r0:b.r0+m-1 and cols b.c0:b.c0+n-1.
  b.r0 = A(1);
  b.m  = A(2);
  b.c0 = A(3);
  b.n  = A(4);
  block = read(fid, 1, char);
  switch(block)
   case 'U' % Outer-product block.
    k   = read(fid, 1,        int );  % Rank.
    b.U = read(fid, [b.m k],  real);  % U stored by column.
    b.V = read(fid, [k b.n],  real)'; % V is stored as V^T: hence V^T is
                                      % stored by col.
    b.B = [];
    nnz = (b.m + b.n)*k;              % Number of nonzeros.
   case 'B' % Dense block.
    b.U = [];
    b.V = [];
    b.B = read(fid, [b.m b.n], real); % Dense block stored by col.
    nnz = b.m*b.n;
   otherwise
    error('Illegal block code.');
  end
  
function [int real char] = FilePrecision(realp)
% realp is 1 (single) or [2] (double)
  int = 'int64';
  if (nargout == 1) return; end
  if (realp == 1)
    real = 'float32';
  else
    real = 'double';
  end
  char = 'uchar';
  
function write(fid,data,type)
  n = fwrite(fid,data,type);
  if (n ~= length(data))
    error(sprintf('Wrote %d instead of %d.',n,length(data)));
  end
  
function A = read(fid,sz,type)
  [A n] = fread(fid,sz,type);
  if (n ~= prod(sz))
    error(sprintf('Read %d instead of %d.',n,prod(sz)));
  end
  
function cfn = RealTypeFn(realp)
  if (realp == 1)
    cfn = @single;
  else
    cfn = @double;
  end

function v7tov8(fnin,fnout)
% Convert a <=v0.7 H-matrix file to a >=v0.8 one.
  ifid = fopen(fnin,'rb');
  int = FilePrecision();
  m  = read(ifid,1,int);
  n  = read(ifid,1,int);
  realp = read(ifid,1,int);
  p  = read(ifid,m,int);
  q  = read(ifid,n,int);
  nb = read(ifid,1,int);
  ofid = fopen(fnout,'wb');
  WriteHeader(ofid,p+1,q+1,realp,nb,-1,-1,-1);
  for (i = 1:nb)
    b = ReadBlock(ifid,realp);
    WriteBlock(ofid,realp,b.r0+1,b.m,b.c0+1,b.n,b.U,b.V',b.B);
  end
  fclose(ifid);
  fclose(ofid);
  
% ------------------------------------------------------------------------------
% Misc

function [nnz n] = ApplyBlockFn(savefn,fn)
  fid = fopen(Filename(savefn),'rb');
  bs = [];
  nnz = 0;
  [p q m n nb realp] = ReadHeader(fid);
  [int real] = FilePrecision(realp);
  for (i = 1:nb)
    [b bnnz] = ReadBlock(fid,realp);
    fn(b);
    nnz = nnz + bnnz;
  end
  fclose(fid);

function fn = FroNorm(savefn)
  global gfn2;
  gfn2 = 0;
  ApplyBlockFn(savefn,@FroNorm2Block);
  fn = sqrt(gfn2);
  clear global gfn2;
  
function FroNorm2Block(b)
  global gfn2;
  if (~isempty(b.B))
    gfn2 = gfn2 + sum(b.B(:).^2);
  else
    gfn2 = gfn2 + sum(sum(b.V'.*((b.U'*b.U)*b.V')));
  end

function [m n realp] = HmatSize(savefn)
  fid = fopen(Filename(savefn),'rb');
  int = FilePrecision();
  m  = read(fid,1,int);
  n  = read(fid,1,int);
  realp = read(fid,1,int);
  fclose(fid);

function in = HmatInfo(savefn)
  fid = fopen(Filename(savefn),'rb');
  int = FilePrecision();
  in.m  = read(fid,1,int);
  in.n  = read(fid,1,int);
  in.realp = read(fid,1,int);
  in.p  = read(fid,in.m,int) + 1;
  in.q  = read(fid,in.n,int) + 1;
  in.nb = read(fid,1,int);
  in.rerr_method = read(fid,1,int); % 1: mrem; 0: brem
  in.tol = read(fid,1,'double');
  est_fronorm = read(fid,1,'double');
  % This is how to construct tol for files from v0.8 to v0.13.
  if (in.rerr_method == 1 && est_fronorm ~= 14.0)
    in.tol = in.tol * est_fronorm;
  end
  fclose(fid);

function f = BlocksPrecon(id,opts)
% Build a preconditioner based on the full blocks of the H-matrix.
  if (nargin < 2) opts = []; end
  if (isfield(opts,'hm_mvp'))
    hm_mvp_fn = opts.hm_mvp;
  else
    hm_mvp_fn = @hm_mvp;
  end
  % I find that a cutoff value of 1e-3 is about right. At least for now it
  % doesn't make sense to make this a user option.
  if (~isfield(opts,'cutoff')) opts.cutoff = 1e-3; end
  keep_going = true;
  f = [];
  while (keep_going && opts.cutoff > eps)
    [I J S] = hm_mvp_fn('bblocks',id,opts.cutoff);
    m = max(I); n = max(J);
    A = sparse(I,J,S,m,n);
    clear I J S;
    setup.type = 'nofill';
    try
      [f.L f.U] = ilu(A,setup);
      keep_going = false;
    catch
      % A has a zero in the diag. Reduce the cutoff.
      opts.cutoff = opts.cutoff/10;
    end
  end
    
function [ce n1 in1] = CondEst1(id,opts)
% This routine copies the basic procedure in Matlab's condest1.
  if (nargin > 1 && isfield(opts,'hm_mvp'))
    hm_mvp_fn = opts.hm_mvp;
  else
    hm_mvp_fn = @hm_mvp;
  end

  n1 = hm_mvp_fn('norm1',id);
  p = BlocksPrecon(id,struct('hm_mvp',hm_mvp_fn));
  tmp = p.L';
  p.L = p.U';
  p.U = tmp;
  clear tmp;
  t = 1;
  in1 = normest1(@(flag,x)CondEst1_InvMvp(flag,x,id,hm_mvp_fn,p),t);
  
  ce = n1*in1;
  %fprintf(1,'hm::CondEst1: n1: %e  in1: %e  ce: %e\n',n1,in1,ce);
  
function f = CondEst1_InvMvp(flag,x,id,hm_mvp,p)
  switch (flag)
   case 'dim'
    f = hm_mvp('getn',id);
   case 'real'
    f = 1;
   case 'notransp'
    [f o] = SolveBxb(id,x,1e-6,...
		     struct('precon',p,'dir_mvp','mvp','hm_mvp',hm_mvp));
    %fprintf(1,'fwd: #mvp = %d\n',o.nmvp);
   case 'transp'
    [f o] = SolveBxb(id,x,1e-6,...
		     struct('precon',p,'dir_mvp','tmvp','hm_mvp',hm_mvp));
    %fprintf(1,'transp: #mvp = %d\n',o.nmvp);
  end

function [ce re] = CondEstFro(id,opts)
  if (nargin < 2) opts = []; end
  if (~isfield(opts,'hm_mvp')) opts.hm_mvp = @hm_mvp; end
  nf2 = opts.hm_mvp('fronorm2',id);
  [nfi re] = NormInvEstFro(id,opts);
  ce = sqrt(nf2)*nfi;
  
function [ne re] = NormInvEstFro(id,opts)
% It's easy to get an estimate of the Frobenius norm of a matrix by performing a
% few MVP with random vectors drawn from certain distributions. One such
% distribution is IID N(0,1) for each element, i.e., randn(n,1). Use this idea,
% combined with SolveBxb, to estimate the Frobenius-norm condition number of
% inv(B).
  if (nargin < 2) opts = []; end
  if (~isfield(opts,'reltol')) opts.reltol = 0.01; end
  if (~isfield(opts,'hm_mvp')) opts.hm_mvp = @hm_mvp; end
  % Each call to SolveBxb of course uses multiple MVP, so we could
  % potentially use quite a number. Still, not a big deal since all it
  % requires is a bunch of H-matrix MVP.
  if (~isfield(opts,'max_nits')) opts.max_nits = opts.hm_mvp('getn',id); end
  if (~isfield(opts,'precon')) opts.precon = 1; end
  
  n = opts.hm_mvp('getn',id);

  if (~isstruct(opts.precon) && opts.precon)
    tt = tic;
    opts.precon = BlocksPrecon(id,opts);
    fprintf(1,'et BlocksPrecon %f\n',toc(tt));
    tmp = opts.precon.L';
    opts.precon.L = opts.precon.U';
    opts.precon.U = tmp;
    clear tmp;
  end
  if (~isfield(opts,'precon') || ~isstruct(opts.precon))
    opts.precon = [];
  end
  
  Ninc = 3;
  N = min(Ninc,opts.max_nits);
  Np = 0;
  nes = [];
  % Loop until the jackknife variance says our estimate is good.
  while (true)
    for (i = 1:(N-Np))
      % MVP of inv(B) with a random vector.
      u = randn(n,1);
      [y o] = SolveBxb(...
	  id,u,1e-6,struct('precon',opts.precon,'hm_mvp',opts.hm_mvp));
      fprintf(1,'%d\n',o.nmvp);
      nes(end+1) = sum(y.^2);
    end
    [ne sigma2] = Jackknife(@(X)mean(X),nes);
    re = sqrt(sigma2)/ne;
    fprintf(1,'%d: %e %f\n',N,sqrt(ne),100*re);
    if (re <= opts.reltol) break; end
    Np = N;
    N = N + Ninc;
    if (N > opts.max_nits) break; end
  end
  ne = sqrt(ne);
  
function [x o] = SolveBxb(id,b,tol,opts)
  if (nargin < 4) opts = []; end
  if (isfield(opts,'hm_mvp'))
    hm_mvp_fn = opts.hm_mvp;
  else
    hm_mvp_fn = @hm_mvp;
    opts.hm_mvp = hm_mvp_fn;
  end
  % Use B(is,is) rather than the full B.
  if (~isfield(opts,'is')) opts.is = []; end
  if (~isfield(opts,'precon')) opts.precon = []; end
  % Forward or transpose problem.
  if (~isfield(opts,'dir_mvp')) opts.dir_mvp = 'mvp'; end
  % Optionally provide a starting point.
  if (~isfield(opts,'x0')) opts.x0 = zeros(length(b),1); end
  if (~isfield(opts,'restart')) opts.restart = 20; end

  opts.n = hm_mvp_fn('getn',id);
  n = opts.n;

  % Save the hm_mvp state associated with the index set for speed.
  if (~isempty(opts.is)) hm_mvp_fn('ststate',id); end
  
  Afun = @(x)SolveBxb_Afun(id,x,opts);

  % Set up and call GMRES
  restart = opts.restart;
  maxit = floor(n/restart);
  if (~isempty(opts.precon))
    if (isa(opts.precon,'struct'))
      if (isfield(opts.precon,'R'))
	precon = @(x)SolveBxb_Mfun_full(opts.precon,x);
      else
	precon = @(x)SolveBxb_Mfun_ILU(opts.precon,x);
      end
    elseif (isa(opts.precon,'function_handle'))
      precon = opts.precon;
    else
      error('opts.precon is not correct.');
    end
  else
    precon = [];
  end
  warning('off','MATLAB:gmres:tooSmallTolerance');
  [x o.flag o.relres o.iter o.resvec] =...
      gmres(Afun,b,restart,tol,maxit,precon,[],opts.x0);
  warning('on','MATLAB:gmres:tooSmallTolerance');
  o.nmvp = restart*o.iter(1) + o.iter(2);
  
  if (~isempty(opts.is)) hm_mvp_fn('strelease',id); end
  
function y = SolveBxb_Afun(id,x,opts)
  if (~isempty(opts.is))
    x1 = zeros(opts.n,1);
    x1(opts.is) = x;
    y = opts.hm_mvp(opts.dir_mvp,id,x1,opts.is,opts.is);
    y = y(opts.is);
  else
    y = opts.hm_mvp(opts.dir_mvp,id,x);
  end
  
function x = SolveBxb_Mfun_full(f,b)
  x = f.Q*(f.U \ (f.L \ (f.P*(f.R \ b))));

function x = SolveBxb_Mfun_ILU(f,b)
  x = f.U \ (f.L \ b);

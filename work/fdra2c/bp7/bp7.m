function varargout = mvd (varargin)
  [varargout{1:nargout}] = feval(varargin{:});
end

function e = env (e)
  e.work_dir = '/home/rverwijs/fdra2c/work';
  e.scratch_dir = '/nfs/camcat001/roos/fdra2c/bp7/sliplaw/ncomp1/dc053/domain_400x400';
  e.suffix='';
  if ~isfield(e,'prefix') e.prefix='ls'; end
end

function addpaths ()
  e = env();
  addpath([e.work_dir '/fdra2c/matlab']);
  addpath([e.work_dir '/dc3dm/matlab']);
%   addpath([e.work_dir])
end

function s = bp7_setup (run_cmd)
% Set up a simulation.
  if (nargin < 1) run_cmd = 1; end
  addpaths()
  vpl=1e-9;
  a0=0.004;
  amax=0.016;
  b=0.01;
  W=400;
  R=200;
  dx=10;
  %wi=12e3;
  %hi=1.5e3;
  widfac=2; %x-width will be l*wid_fac  
  T = 1; %duration of nucleation stress perturbation 
  Rnuc = 150;
  tau0 = 1.75e6;
  xn = -50; % hypocenter for nucleation stress pert
  yn = -50;

  o = bp7_setopts(a0, amax, b, W, R, vpl, T, Rnuc, tau0, xn, yn);
  o.ncomp = 1; %TODO set to 2. At the moment this won't load BC.
  o.dx = dx;
  o.widfac = widfac;
  addpaths(); 
  e=env();
  sourceFile = mfilename('fullpath');
%   disp([sourceFile '.m'])
  status = copyfile([sourceFile '.m'], e.scratch_dir);
  copyfile /home/rverwijs/fdra2c/work/fdra2c/bp7/readme e.scratch_dir 
  copyfile /home/rverwijs/fdra2c/work/fdra2c/bp7/script.sh e.scratch_dir
  e.description='SEAS BP7 as described in the pdf from 4/2023.';

  fout = fopen([e.scratch_dir '/readme'],'a');
  fprintf(fout, '%s', e.description);
  fclose(fout);
  
  if (nargin>1 & ~isempty(init_cond)) o.init_cond=init_cond; end
 
  scriptfile = fopen([e.scratch_dir,'/','script.sh'],'a');

  o.tf = 100;
  o.vinit=1e-9;
  o.vinit_nuc=1e-3;
  o.vzero=1e-20;
  o.sv=[1 0 0]; %slip boundary condition vector.

  o.ffn = [o.dir '/bp7_slip_1t' num2str(o.ncomp) 'd_' num2str(dx)  'm'];
  o.hmfn = [o.dir '/bp7_slip_1t' num2str(o.ncomp) 'd_' num2str(dx) 'm_hm'];

  iw=1; %I'm not sure if s.hm needs to be a cell so I'll keep it this way.

  s.hm{iw}.cb = bp7_write_build_kvf(o);
  s.hm{iw}.cm = bp7_write_mesh_kvf_uni(o); %this produces uniform mesh.
  if (run_cmd)
    mysystem(s.hm{iw}.cm.cmd);
    mysystem(s.hm{iw}.cb.cmd);
  else
    fprintf(scriptfile, s.hm{iw}.cm.cmd);
    fprintf(scriptfile, s.hm{iw}.cb.cmd);
  end

  %Four structures are greens functions between elements, two are from BCs. Note: If "compress" takes really long it may be best to pass the slip vector orientation to fdra and do the calculation after, since the green's function to to any slip BC is a linear combination of the ones for slip in the [1,0,0] and [0,1,0] directions.
  s.hm{iw}.cc{1,1} = bp7_write_compress_kvf(o,1,1);
  s.hm{iw}.bc{1} = bp7_write_compress_kvf_bc(o,o.sv,1);

  if o.ncomp==2
    s.hm{iw}.cc{1,2} = bp7_write_compress_kvf(o,1,2);
    s.hm{iw}.cc{2,1} = bp7_write_compress_kvf(o,2,1);
    s.hm{iw}.cc{2,2} = bp7_write_compress_kvf(o,2,2);
    s.hm{iw}.bc{2} = bp7_write_compress_kvf_bc(o,o.sv,2);
  end

  for i=1:o.ncomp
  if (run_cmd)
       mysystem(s.hm{iw}.bc{i}.cmd);
       system(['rm ' s.hm{iw}.bc{i}.hm_write_filename '.hm']);      

       if o.ncomp==1 system(['ln -s ' o.hmfn '_comp1bc.bc ' o.hmfn '.bc']); 
       end
  else
       fprintf(scriptfile, s.hm{iw}.bc{i}.cmd);
    end

    for j=1:o.ncomp
       if (run_cmd)
	    mysystem(s.hm{iw}.cc{i,j}.cmd);
	    system(['rm ' s.hm{iw}.cc{i,j}.hm_write_filename '.bc'])
       else
            fprintf(scriptfile, s.hm{iw}.cc{i,j}.cmd);
       end
    end
  end
%  try
   for (im = 1)
    s.fs{iw, im}.cf = bp7_write_fdra_kvf(o);
%   end

%   s.hm{iw}.cm = bp7_write_mesh_kvf(o, s.fs{iw, im}.cf);

   fprintf(scriptfile, s.fs{iw}.cf.cmd);3
  end

  fclose(scriptfile);

  if (run_cmd)
    system(['ln -s ' o.hmfn '_comp11.hm ' o.hmfn '_comp11.hmat']);
    system(['ln -s ' o.hmfn '_comp12.hm ' o.hmfn '_comp12.hmat']);
    system(['ln -s ' o.hmfn '_comp21.hm ' o.hmfn '_comp21.hmat']);
    system(['ln -s ' o.hmfn '_comp22.hm ' o.hmfn '_comp22.hmat']);
  end
  
  fname=[o.ffn '.mat'];
  save(fname,'s');
end

function o = bp7_setopts (a0, amax, b, W, R, vpl, T, Rnuc, tau0, xn, yn)
  s2y = 1/(365.25*24*60*60);
  e = env();
  o.dir = e.scratch_dir;
  o.prefix = e.prefix;
  o.tol = 1e-6;
  % Number of CPU cores to use:
  o.nthreads = 1; 
  o.output_factor = 10;
  o.rel_tol = 1e-3;
  
  o.evolution = 'slip';

  rho = 2670;
  vs = 3.464e3;
  
  o.mu = vs^2*rho;
  o.nu = 0.25;
  o.vs = vs;
  o.v_creep = vpl;
  o.sigma = 25e6;
  o.d_c = 0.53e-3;
  o.a0 = a0;
  o.amax = amax;
  o.b = b;
  o.W = W;
  o.R = R;
  o.T = T;
  o.xn = xn;
  o.yn = yn;
  o.Rnuc = Rnuc;
  o.tau0 = tau0;

%  o.l = l;
%  o.H = H;
%  o.wi = wi;
%  o.hi = hi;
  
  o.models = {'neutral'};
  o.model = o.models{1};
end

function bp7_write_script (s, script_fn)
  fid = fopen(script_fn, 'w');
  for (iw = 1:numel(s.hm))
    fprintf(fid, '%s', s.hm{iw}.cc.cmd);
    for (im = 1 : size(s.fs, 2))
      fprintf(fid, '%s', s.fs{iw, im}.cf.cmd);
    end
    fprintf(fid, '\n');
  end
  fclose(fid);
end

function c = bp7_write_mesh_kvf_uni (o)
 % The mesh for dc3dm.
   c.mesh_write_filename = o.hmfn;
   c.max_len = o.dx;
   c.y = o.W*[-1 1];
   c.x = o.W*[-1 1];
   c.do_uniform = 1;
   c = write_kvf_finish(o.hmfn, c, 'mesh');
end

function c = bp7_write_build_kvf (o)
% Boundary conditions, half- or fullspace, and other parameters to make the
% dislocation-traction operator.
  bfn = o.hmfn;
  c.mesh_read_filename = bfn;
  c.build_write_filename = bfn;
  c.do_fullspace = 1;
  c.neighborhood = 0;
 
  %east and west: periodic BC. Not dominant (but doesn't really matter for periodic).
%  c.ewpbc = 1;
%  c.nvbc = 2;
%  c.svbc = 3;

  c.evbc = 1;
  c.wvbc = 4;
  c.nvbc = 2;
  c.svbc = 3;
  c = write_kvf_finish(o.hmfn, c, 'build');
end

function c = bp7_write_compress_kvf (o,sourcedir,recdir)
% Final parameters for the compressed operator.
  bfn = [o.hmfn '_comp' num2str(sourcedir), num2str(recdir)];
  c.mu = o.mu; c.nu = o.nu;
  c.build_read_filename = o.hmfn;
  c.tol = o.tol;
  c.allow_overwrite = 1;
  c.nthreads = o.nthreads;

  %Write out files for required component:
  c.src_disl = [0 0 0];
  c.rcv_traction = c.src_disl;
  c.src_disl(sourcedir)=1;
  c.rcv_traction(recdir)=1;

  c.hm_write_filename = bfn;
  c = write_kvf_finish(bfn, c, 'compress');

end

function c = bp7_write_compress_kvf_bc (o,sourcevector,recdir)
% Final parameters for the compressed operator. This is used for boundary conditions.
  bfn = [o.hmfn '_comp' num2str(recdir) 'bc'];
  c.mu = o.mu; c.nu = o.nu;
  c.build_read_filename = o.hmfn;
  c.tol = o.tol;
  c.allow_overwrite = 1;
  c.nthreads = o.nthreads;

  %Write out files for required component:
  c.src_disl = sourcevector;
  c.rcv_traction = [0 0 0];
  c.rcv_traction(recdir)=1;

  c.hm_write_filename = bfn;
  c = write_kvf_finish(bfn, c, 'compress')

end

function c = bp7_write_fdra_kvf (o)
  e = env();
  s2y = 1/(365.25*24*60*60);

  c = o;
  c.rmesh_filename = o.hmfn;
  c.rs = dc3dm.mRects(c.rmesh_filename);
  [c.xc c.yc] = dc3dm.mCC(c.rs);
  c.nelem = numel(c.xc);

  one = ones(size(c.xc));
  c.v0 = 1e-6;
  mu0=0.6;
  c.mu0 = mu0*one;
  c.b = o.b*one;
  c.d_c = o.d_c*one;
  c.s_normal = o.sigma*one;

  % Set a values:
  r = sqrt(c.yc.^2+c.xc.^2);
  out = r > o.R;
  in = r <= o.R;

  c.a(in) = o.a0;
  c.a(out) = o.amax;

  c.a=c.a(:);

  c.evolution = o.evolution;
  c.use_vcutoff = 0;
  c.ncomp = o.ncomp;
  
  v_s = o.vs;
  c.eta = 0.5*o.mu/v_s;
  
  c.stress_fn = 'h-matrix-nucl';
  %c.stress_fn = 'h-matrix';
  c.v_creep = o.v_creep;

  % Set nucleation variables:
  % TODO shift by y2, y3
  c.T = o.T;
  c.tau0 = o.tau0;
  rn = sqrt((c.yc-c.yn).^2+(c.xc-c.xn).^2);
  c.G1 = exp(rn.^2./(rn.^2-o.Rnuc^2)) .* (rn<o.Rnuc);

  c.ti = 0;
  c.tf = o.tf/s2y();
  
  c.hm_filename = o.hmfn;
  c.rel_tol = o.rel_tol;
  fac = o.output_factor;
  c.save_filename = o.ffn;
  c.disp_every = 10;
  c.stop_indicator = 'stop.ind';
  c.stop_check_frequency = 100*fac;
  c.save_v_every = fac;
  c.save_slip_every = fac;
  c.save_state_every = fac;
  c.allow_overwrite = 1;
 
  %Region of higher prestress and slip rates:
%  nxc = -o.l/2+o.hi+o.wi/2;
%  nyc = o.H-o.hi-o.wi/2;

%  nucreg = c.xc>nxc-o.wi/2 & c.xc<nxc+o.wi/2  & c.yc>nyc-o.wi/2 & c.yc<nyc+o.wi/2;

  c.v_init = [o.vinit*one];%; o.vzero*one];
%  c.v_init(nucreg)=o.vinit_nuc;

  % define tau00 = tau0-etaV
  tau00 = c.sigma*c.a.*asinh(c.v_init/2/c.v0.*exp((mu0+o.b*log(c.v0/o.vinit))./c.a));

  gamma_init = c.a./c.b .* log(2*c.v0./c.v_init .* sinh(tau00./c.a/o.sigma)) - mu0/o.b;
  c.chi_init = gamma_init + log(c.v_init/c.v0);

  if o.ncomp==2
      c.v_init = [o.vinit*one; o.vzero*one];
%      c.v_init(nucreg)=o.vinit_nuc;
  end

%  anuc = c.a(nucreg);
%  gamma_init = anuc/o.b .* log(2*c.v0/o.vinit * sinh(1.1*tau00./anuc/o.sigma)) - mu0/o.b;
%  c.chi_init(nucreg) = gamma_init + log(o.vinit/c.v0);

  c.kvf = [o.ffn '_f'];
  c.cmd = sprintf( ...
      ['%s/fdra2c/bin/add_bc_to_fdra_kvf %s.kvf;\n',...
       '%s/fdra2c/bin/fdra %s.kvf;\n'], ...
       e.work_dir, c.kvf, e.work_dir, c.kvf);
       %       'time mpirun -np %d %s/fdra2c/bin/fdra %s.kvf;\n'],...
       %   e.work_dir, c.kvf, o.nthreads, e.work_dir, c.kvf);
  dc3dm.WriteKvf(c.kvf, c, true);
end

% ------------------------------------------------------------------------------
% Utils.

function sd = transfer_fields (sd, ss, flds)
  if (~iscell(flds)) flds = {flds}; end
  for (i = 1:numel(flds))
    if (isfield(ss, flds{i})) sd.(flds{i}) = ss.(flds{i}); end
  end
end

function mysystem (varargin)
  cmd = sprintf(varargin{:});
  %pr([cmd '\n']);
  %[s r] =...
      system(cmd);
end

function ck = write_kvf_finish (bfn, ck, command)
  e = env();
  ck.command = command;
  ck.kvf = [bfn '_' command(1)];
  ck.cmd = sprintf('time %s/dc3dm/bin/dc3dm %s.kvf;\n', e.work_dir, ck.kvf);
  dc3dm.WriteKvf(ck.kvf, ck, 1);
end

function axisall ()
  axis equal; axis xy; axis tight; 
end
function o = dopt (o, fld, val)
  if (~isfield(o, fld)) o.(fld) = val; end
end
function s = rmsp (s)
  s(s == ' ') = [];
end

function hs = calc_hstar_bma (mu, nu, a, b, d_c, sigma)
  hs = pi*mu/(1 - nu)*d_c ./ (4*sigma.*(b - a));
end
function hs = calc_hstar_b (mu, nu, b, d_c, sigma)
  hs = 1.377*mu/(1 - nu)*d_c ./ (sigma.*b);
end

function varargout = mvd (varargin)
  [varargout{1:nargout}] = feval(varargin{:});
end

function e = env ()
  e.work_dir = '/home/ambrad/basic3d/work';
  e.scratch_dir = '/scratch/ambrad/b3d';
end

function addpaths ()
  e = env();
  addpath([e.work_dir '/fdra2c/matlab']);
  addpath([e.work_dir '/dc3dm/matlab']);
end

function s = cl_setup (run_cmd)
% Set up a simulation.
  if (nargin < 1) run_cmd = 1; end
  o = cl_setopts();
  % Crack radius. For demo purposes, look at just a very small crack:
  rads = [1];
  for (iw = 1:numel(rads))
    o.asp_radius = rads(iw);
    o.dmn_width = 2*rads(iw) + 10;
    o.output_factor = rads(iw)*3;
    o.tf = rads(iw)/2;
    o = cl_update_opts(o);
    s.hm{iw}.cm = cl_write_mesh_kvf(o);
    s.hm{iw}.cb = cl_write_build_kvf(o);
    if (run_cmd)
      mysystem(s.hm{iw}.cm.cmd);
      mysystem(s.hm{iw}.cb.cmd);
    end
    s.hm{iw}.cc = cl_write_compress_kvf(o);
    for (im = 1)
      s.fs{iw, im}.cf = cl_write_fdra_kvf(o);
    end
  end
end

function o = cl_setopts ()
  e = env();
  o.dir = e.scratch_dir;
  o.rfac = 1;
  o.tol = 1e-6;
  % Number of CPU cores to use:
  o.nthreads = 8; % Use 16 if on redoubt, doing serious work, and nobody else is
                  % using the machine; 8 for this example because I saw another
                  % user on redoubt when writing it.
  o.output_factor = 10;
  o.rel_tol = 1e-3;
  o.tf = 1;
  
  o.evolution = 'aging';
  o.mu = 3e10;
  o.nu = 0;
  o.v_creep = 1e-9;
  % Nothing subtle seems to be going on: if I increase this factor to f, then
  % events happen ~f times more frequently.
  o.v_creep_fac = 1;
  o.sigma = 50e6;
  o.d_c = 1e-4;
  o.a = 0.015;
  % Setting to o.a makes a/b = 0.5.
  o.ab_delta = o.a; 0.004;
  o.asp_radius = 2; % [h*_{b-a}]
  o.dmn_width = 5;
  
  o.models = {'neutral'};
  o.model = o.models{1};
end

function cl_write_script (s, script_fn)
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
function o = cl_update_opts (o)
  % cl_ was the original with a/b = 0.5. h*_{b-a} was wrong.
  % cl2_ was for a/b used in C&L.
  % cl3_ is back to a/b = 0.5 to study scaling at large w.
  o.ffn = sprintf( ...
    rmsp('%s/cl3_ m%s ar%1.1f dw%d _rf%d'), ...
    o.dir, o.model(1:2), o.asp_radius, o.dmn_width, o.rfac);
  o.hmfn = sprintf( ...
    rmsp('%s/cl3_ ar%1.1f dw%d _rf%d tol%d'), ...
    o.dir, o.asp_radius, o.dmn_width, o.rfac, ...
    round(log10(o.tol)));
end

function c = cl_write_mesh_kvf (o)
% The mesh for dc3dm.
  c.mesh_write_filename = o.hmfn;
  c.max_len = calc_hstar_b(o.mu, o.nu, o.a + o.ab_delta, o.d_c, o.sigma)/ ...
      (5*o.rfac);
  h = calc_hstar_bma(o.mu, o.nu, o.a, o.a + o.ab_delta, o.d_c, o.sigma);
  c.y = 0.5*o.dmn_width*h*[-1 1];
  c.x = c.y;
  c.do_uniform = 1;
  c = write_kvf_finish(o.hmfn, c, 'mesh');
end
function c = cl_write_build_kvf (o)
% Boundary conditions, half- or fullspace, and other parameters to make the
% dislocation-traction operator.
  bfn = o.hmfn;
  c.mesh_read_filename = bfn;
  c.build_write_filename = bfn;
  c.do_fullspace = 1;
  c.neighborhood = 0;
  dir = 'enws';
  for (i = 1:4) c.([dir(i) 'vbc']) = i; end
  c = write_kvf_finish(o.hmfn, c, 'build');
end
function c = cl_write_compress_kvf (o)
% Final parameters for the compressed operator.
  bfn = o.hmfn;
  c.src_disl = [1 0 0];
  c.rcv_traction = c.src_disl;
  c.mu = o.mu; c.nu = o.nu;
  c.build_read_filename = bfn;
  c.tol = o.tol;
  c.hm_write_filename = bfn;
  c.allow_overwrite = 0;
  c.nthreads = o.nthreads;
  c = write_kvf_finish(o.hmfn, c, 'compress');
end

function c = cl_write_fdra_kvf (o)
  e = env();

  if (~exist([o.hmfn '_comp11.hmat'], 'file'))
    mysystem(sprintf('ln -s %s.hm %s_comp11.hmat;', o.hmfn, o.hmfn));
  end

  c = o;
  c.rmesh_filename = o.hmfn;
  c.rs = dc3dm.mRects(c.rmesh_filename);
  [c.xc c.yc] = dc3dm.mCC(c.rs);
  c.nelem = numel(c.xc);

  one = ones(size(c.xc));
  c.v0 = 1e-6;
  c.mu0 = 0.6*one;
  c.a = o.a*one;
  c.b = (o.a + o.ab_delta)*one;
  c.d_c = o.d_c*one;
  c.s_normal = o.sigma*one;
  
  r = sqrt(c.xc.^2 + c.yc.^2);
  h = calc_hstar_bma(o.mu, o.nu, o.a, o.a + o.ab_delta, o.d_c, o.sigma);
  mask = r > o.asp_radius*h;
  c.b(mask) = o.a;
  c.a(mask) = o.a + o.ab_delta;
  
  c.evolution = o.evolution;
  c.use_vcutoff = 0;
  c.ncomp = 1;
  
  v_s = 3.7e3;
  c.eta = 0.5*o.mu/(1 - o.nu)/v_s;
  
  c.stress_fn = 'h_matrix';
  c.v_creep = o.v_creep_fac*o.v_creep;

  c.ti = 0;
  c.tf = o.tf/s2y();
  
  c.hm_filename = o.hmfn;
  c.rel_tol = o.rel_tol;
  fac = o.output_factor;
  c.save_filename = o.ffn;
  c.disp_every = 1000;
  c.stop_indicator = 'stop.ind';
  c.stop_check_frequency = 100*fac;
  c.save_v_every = fac;
  c.save_slip_every = fac;
  c.save_state_every = 10*fac;
  c.allow_overwrite = 1;
  
  c.v_init = c.v_creep*one;
  c.chi_init = zeros(size(c.a));
    
  c.kvf = [o.ffn '_f'];
  c.cmd = sprintf( ...
      ['%s/fdra2c/bin/add_bc_to_fdra_kvf %s.kvf;\n',...
       'time mpirun -np %d %s/fdra2c/bin/fdra %s.kvf;\n'],...
    e.work_dir, c.kvf, o.nthreads, e.work_dir, c.kvf);
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

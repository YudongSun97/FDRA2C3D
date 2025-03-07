function varargout = mvd (varargin)
% function s = q2fdra2c (fn, run_hmcomp)
% converts between matlab and C++ FDRA.
% copies parameters and initial conditions from structures used by the matlab version of FDRA into structures and key-value files (KVF) used by the C++ version of FDRA.
% also reads full Green's functions matrix from original simulation and performs H-matrix compression.
% 
% fn = name of original simulation (matlab FDRA) 
% run_cmd = flag indicating if H-matrix compression should be performed at run time. This uses a system call. if run_cmd=0 you will need to run the command on your own before calling FDRA.% s = structure with parameters for C++ simulation (similar to input KVF).

%  [varargout{1:nargout}] = feval(varargin{:});
   varargout{1} = simul_setup(varargin{:});
end

% Set up the environment
function e = env(e)
  % fdra1_dir: directory containing fdra1.33 (matlab FDRA)
  % work_dir: directory containing fdra2c
  % scratch_dir: directory for output files
  % hmmvp_dir: directory containing hmmvp (H-matrix compression)
  e.fdra1_dir = '/home/camcat/code/matlab/fdra1.33-master'; 
  e.work_dir = '/home/camcat/code/fdra2c-extern_gf/work';
  e.scratch_dir = '/nfs/camcat001/camcat/fdra2c/tests';
  e.hmmvp_dir = '/home/camcat/code/fdra2c-extern_gf/work/fdra2c/external/hmmvp';

  e.dec='v2';
  e.suffix='';
  if ~isfield(e,'prefix') e.prefix=''; end
end

function addpaths ()
  e = env();
  addpath([e.work_dir '/fdra2c/matlab']);
  addpath(e.fdra1_dir);
  addpath([e.fdra1_dir '/extras']);
end

% ----------------------------------------------------- %
%                   Set up a simulation                 %
% ----------------------------------------------------- %
function s = simul_setup (fn, run_cmd)
  % if run_cmd=1, run the commands to compress the Green's functions (calls hmmvp)
  if (nargin < 2) run_cmd = 1; end

  % Set up the environment and add paths
  e=env();
  addpaths();
 
  % Load structures from matlab version of fdra
  warning('off', 'all');
  [c,l]=Util('ReadConsts',fn,1,1);
  warning('on', 'all');

  % set up options
  o = simul_setopts(c,l);
 
  % filenames: simulation (ffn) and H-matrices (hmfn)
  [~,i]=find(fn=='/',1,'last');
  fn2 = fn(i+1:end);
  o.save_filename = [e.scratch_dir '/' o.prefix fn2 o.dec];
  o.hm_filename = [e.scratch_dir '/' o.prefix fn2 o.dec '_hm'];

  % load Greens functions matrices
  % The Matlab version of FDRA uses one of two variables to save Green's functions (c.tgf or l.Gs)
  % Ideally this should be standardized
  try 
     s.hm.cc{1,1} = c.tgf;
     if (c.use_Gn) s.hm.cc{1,2} = c.sgf; end
  catch 
     s.hm.cc{1,1} = l.Gs;
     if (c.use_Gn) s.hm.cc{1,2} = l.Gn; end
  end

  % calculate H-matrices from Greens functions matrices
  s.hm.cc{1,1} = gf2hm(s.hm.cc{1,1}, [o.hm_filename '_comp11.hm'], ...
	o.xc, o.yc, 0*o.xc, run_cmd, o.hmmvp_mode);
  if (c.use_Gn) 
        s.hm.cc{1,2} = gf2hm(s.hm.cc{1,2}, [o.hm_filename '_comp1n.hm'], ...
	   o.xc, o.yc, 0*o.xc, run_cmd, o.hmmvp_mode);
  end

  % -------- run commands, write output files ----------- %

  % hard link H-matrices with .hmat extension (compatible with fdra2c)
  % compXY indicates source and receiver components (1=along strike, 2=along dip, n=normal).
  % we're assuming no opening so there is no normal source.
  if (run_cmd)
    system(['ln -s ' o.hm_filename '_comp11.hm ' o.hm_filename '_comp11.hmat']);
    if c.use_Gn system(['ln -s ' o.hm_filename '_comp1n.hm ' o.hm_filename '_comp1n.hmat']); end
    if o.ncomp==2
      system(['ln -s ' o.hm_filename '_comp12.hm ' o.hm_filename '_comp12.hmat']);
      system(['ln -s ' o.hm_filename '_comp21.hm ' o.hm_filename '_comp21.hmat']);
      system(['ln -s ' o.hm_filename '_comp22.hm ' o.hm_filename '_comp22.hmat']);
      if c.use_Gn system(['ln -s ' o.hm_filename '_comp1n.hm ' o.hm_filename '_comp2n.hmat']); end
    end
  end

  % write out readme file
  e.description=['Copy of: ' fn];
  fout = fopen([e.scratch_dir '/readme'],'a');
  fprintf(fout, '%s\n', e.description);  
  fclose(fout);

  % write key-value file and script file
  s.fs.cf = simul_write_fdra_kvf(o);
  scriptfile = fopen([e.scratch_dir,'/','script.sh'],'a');
  if ~run_cmd
    for n=1:length(s.hm.cc)
      fprintf(scriptfile, [s.hm.cc{n}.cmd '\n']);
    end
  end
  fprintf(scriptfile, s.fs.cf.cmd);
  fclose(scriptfile);

  %save structure to matlab file
  fname=[o.save_filename '.mat'];
  save(fname,'s');

end

% ----------------------------------------------------- %
%                   Auxiliary functions                 %
% ----------------------------------------------------- %

function o = set_defaults();

  e=env();
  o.ncomp = 1; %number of slip components. 1 for 1D, 1 or 2 for 2D.
  o.use_Gn = 0; %include greens functions for normal stresses? (not yet implemented in fdra2c)
  o.dec=[e.dec 'no_Gn']; %add suffix to filenames to indicate whether normal stresses are included.
  o.prefix=e.prefix;
  o.dir = e.scratch_dir;

  if o.use_Gn error('Normal stresses not yet implemented in fdra2c'); end

  o.hmmvp_mode = 'omp'; %parallelization mode for hmmvp.
  o.stress_fn = 'h-matrix';
  
  o.tol = 1e-6; % tolerance
  o.rel_tol = 1e-3;
  o.nthreads = 8; % Number of CPU cores to use
  
  output_factor = 10; % How often to save output files
  o.disp_every = 100;
  o.stop_indicator = 'stop.ind';
  o.stop_check_frequency = 100*output_factor;
  o.save_v_every = output_factor;
  o.save_slip_every = output_factor;
  o.save_state_every = output_factor;
  o.allow_overwrite = 1;

end


% ------------- Set options for fdra2c based on structures c and l from matlab fdra --------------- %
function o = simul_setopts (c,l);
  cc = @(x)(x(2:end)+x(1:end-1))/2;
 
  % Set defaults
  o = set_defaults();

  % Extract coordinate arrays
  o.xc = cc(c.x);
  try o.yc = cc(c.y);  %rough fault
  catch o.yc = 0*o.xc; %smooth fault (yc not defined)
  end

  % Copy settings and parameters
  if c.evolution==1 o.evolution = 'aging';
  else o.evolution = 'slip';
  end  

  % Copy elastic parameters
  o.mu = c.G;
  o.nu = c.poisson;
  o.vs = c.Vs;
  o.eta = 0.5*o.mu/o.vs;

  if (o.ncomp==2) 
    o.vzero = 1e-20; %slip direction in secondary direction
  end 

  % Copy friction parameters, normal stress
  o.nelem = numel(o.xc);
  one = ones(size(o.xc));
  o.v0 = c.V0;
  o.mu0 = c.mu_0*one;
  o.d_c = c.Dc;
  o.a = c.a;
  o.b = c.b;
  if isfield(c,'kes') o.kes=c.kes; else o.kes=0; end %minimum normal stress (for rough fault)
  o.s_normal = c.s_normal;

  % Copy flags for additional frictional laws
  o.use_vcutoff = c.use_vcutoff;
  o.use_dynweak = c.flash_heating;
  if (o.use_dynweak)
    o.dw_fw = c.flash_fw*one;
    o.dw_vw = c.flash_vw*one;
    o.dw_n = c.flash_n;
  end  

  % Set loading type.
  % TODO this is quite approximate and it may fail for unexpected loading functions
  if ~isfield(l,'delta_tau_fn') || isempty(l.delta_tau_fn);
     o.load = 'bc'; %boundary conditions
  elseif contains(func2str(l.delta_tau_fn), 'loading')
     o.load = 'uniload';  %uniform loading
  elseif contains(func2str(l.delta_tau_fn), 'resolve')
     o.load = 'resolve';  %resolve a remote stress field %TODO this is not yet implemented in fdra2c
     warning('''Resolve'' loading not yet implemented in fdra2c'); %note that below I use an approximation for this case
  else
     error('Could not interpret loading.')
  end

  % Set parameters associated with loading
  switch o.load
  case 'uniload'
      o.uniload = 1;
      o.taudot = c.taudot;
  case 'resolve' %FIXME implement this properly in fdra2c
      o.uniload = 1;
      o.taudot = c.Sij.sigmaD/2; 
      warning('Approximating resolves stresses by uniform taudot');
  case 'bc'
      error('BC loading has not yet been implemented (fix this!)');
  end

  % fdra2c calculates uniform loading rate as: taudot = taudot_o_vcreep*v_creep
  % where v_creep is the creep velocity and taudot_o_vcreep is the loading rate normalized by the creep velocity.
  if o.uniload
     o.v_creep = 1.0;
     if length(o.taudot)==1 o.tdot_o_vcreep=(o.taudot/o.v_creep)*one;
     else o.tdot_o_vcreep=(o.taudot(:)/o.v_creep);
     end
  else 
     o.v_creep = c.Vpl_creep;
  end

  % start and end time
  o.ti = c.ts;
  o.tf = c.tend;
  
  o.v_init = c.V0*exp(c.psi_init);
  o.chi_init = c.chi_init;
  if o.ncomp==2
      o.v_init = [c.v_init; o.vzero*one];
  end

  % Add normal stress perturbation corresponding to roughness
  if ~isfield(c,'is_rough') c.is_rough=0; end
  if c.is_rough==2 %pseudoroughness using normal stresses as proxy for geometry
    deltasig = y_to_sigma(o.xc, o.yc, 1, o.mu/(1-o.nu)).*c.slip_init';
    o.s_normal = o.s_normal + deltasig(:)';
    o.s_normal = max(o.s_normal, o.kes);
  elseif c.is_rough==1 %real roughness
    error('Roughness not yet implemented in fdra2c');
  end

end


% -------- Write key-value file for fdra2c ----------- %
function cf = simul_write_fdra_kvf (o)
  e=env();
  cf=o;
  cf.kvf = [o.save_filename '_f'];

  if strcmp(o.load,'bc')
      cf.cmd = sprintf( ...
          ['%s/fdra2c/bin/add_bc_to_fdra_kvf %s.kvf;\n',...
          '%s/fdra2c/bin/fdra %s.kvf;\n'], ...
          e.work_dir, cf.kvf, e.work_dir, cf.kvf);
  else
    cf.cmd = sprintf( ...
          ['%s/fdra2c/bin/fdra %s.kvf;\n'], e.work_dir, cf.kvf);
  end
  kvf('Write',[cf.kvf '.kvf'], o, true);
end

% convert Green's function matrix to structure used by hmmvp
function c = gf2hm(G, fname, xc, yc, zc, run, mode);
% function c = gf2hm(G, fname, xc, yc, zc, run, mode);
% G = Green's function matrix
% fname = filename for compressed matrix and kvf file
% xc, yc, zc: source points (assumed same as receivers)
% run: flag indicating if the system command should be executed. 
% mode: hmmvp mode (can be 's', 'omp', 'mpi' for serial vs parallel).
%
% c is a structure used by hmmvp to compress the matrix.
% See hmmvp documentation for description of required fields (more specifically, type: ./bin/hmmvpbuild_<mode> help compress where <mode>= indicates parallelization.

   if ~exist('run','var') run=0; end
   if ~exist('mode','var') mode='omp'; end
   
   c.eta=3;
   c.err_method='brem-fro';
   %Tolerance. These are the values from FDRA:
   %[-6, -9, -12, -16]={very loose, loose, moderate, tight}.
   c.tol = 1e-6;
   c.nthreads=8;
   c.command='compress';
   c.allow_overwrite=1;
   c.greens_fn='external';
   
   c.write_hmat_filename=fname;
   c.kvf=[c.write_hmat_filename '.kvf'];
   c.G=G(:,1:end-1);
   c.Nr=length(xc);
   c.Ns=length(xc);
   c.xr=[xc(:), yc(:), zc(:)]';
   c.xs=[xc(:), yc(:), zc(:)]';
   
   e=env();
   hmmvp_bin=[e.hmmvp_dir '/bin/hmmvpbuild_' mode];
   %Save command to compress matrices:
   c.cmd = [hmmvp_bin ' ' c.kvf];

   kvf('Write', c.kvf, c, 1);
   if (run) 
    disp(c.cmd);
    system(c.cmd); 
  end

end

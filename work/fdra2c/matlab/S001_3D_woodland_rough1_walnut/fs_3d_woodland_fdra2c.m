function varargout = mvd (varargin)
    % function s = fs_2d_multifault_q2fdra2c (run_hmcomp)
    % run_cmd = flag indicating if H-matrix compression should be performed at run time. This uses a system call. if run_cmd=0 you will need to run the command on your own before calling FDRA.% s = structure with parameters for C++ simulation (similar to input KVF).
    % Yudong: aka run_hmcomp
    %  [varargout{1:nargout}] = feval(varargin{:});
       varargout{1} = simul_setup(varargin{:});
end

% Set up the environment
function e = env(e)
    % work_dir: directory containing fdra2c
    % scratch_dir: directory for output files
    % hmmvp_dir: directory containing hmmvp (H-matrix compression)
  
    e.work_dir = '/home/users/yudongs/FDRA2C/fdra2c/work';
    e.scratch_dir = '/home/users/yudongs/fdra2c/tests001';
    e.hmmvp_dir = '/home/users/yudongs/FDRA2C/fdra2c/work/fdra2c/external/hmmvp';
    
    e.dec='v2';
    e.suffix='';
    if ~isfield(e,'prefix') e.prefix=''; end
end

function s = simul_setup (run_cmd)
    % Yudong: 2025-06-23
    % This function is used to setup the simulation environment for the 2D multifault simulation.
    % if run_cmd=1, run the commands to compress the Green's functions (calls hmmvp)
    if (nargin < 1) run_cmd = 1; end

    % Set up the environment and add paths
    e=env();
    addpaths();

    % Yudong: 2026-06-12
    % load geometry and Green's functions.
    load(fullfile(e.scratch_dir, 'wl_rough1_walnut.mat'));

    cc = @(x)(x(2:end)+x(1:end-1))/2;
 
    % Set defaults
    o = set_defaults();

    % Yudong: 2025-06-25
    % just write the parameters

    % flags 
    o.evolution = 'aging'; % 'aging' or 'slip'
    o.use_Gn = 1; %include greens functions for normal stresses
    if o.use_Gn ==1
        o.dec=[e.dec 'Gn']; %add suffix to filenames to indicate whether normal stresses are included.
    end
    fn2 = ''; % the name of the simulation
    o.save_filename = [e.scratch_dir '/' o.prefix fn2 o.dec];
    o.hm_filename = [e.scratch_dir '/' o.prefix fn2 o.dec '_hm'];
    o.use_vcutoff = 0; %default: do not use velocity cutoff
    o.use_dynweak = 0; %default: do not use dynamic weakening (flash heating)
    if (o.use_dynweak)
        error('Dynamic weakening (flash heating) is not yet implemented in fdra2c');
    end
    if (o.ncomp==2) % if consider the slip in the secondary direction
        o.vzero = 1e-20; %slip direction in secondary direction
    end 
    
    % elastic parameters in the bulk
    o.mu = 3.0e10; % Shear modulus [Pa]
    o.nu = 0.25; % Poisson ratio    
    o.vs = 3.7e3; % S wave velocity [m/s]
    o.eta = 0.5*o.mu/o.vs;  % radiation damping parameter
    Gprime = o.mu/(1 - o.nu);   % G/(1-nu) for plane strain; G for anti-plane
    
    % friction parameters on the interface
    a = 0.011;
    b = 0.02;
    Dc = 0.968e-6/2;

    % other parameters
    o.v0 = 1e-8;                % Reference velocity [m/s]
    Vpl_creep = 0;      % set to 0 when only using uniform or resolved loadings
    o.mu_0 = 0.6; % friction coefficient
    sigma = 10e6;   % normal stress [Pa]


    % lengths & spacing 
    L_over_Ln = 10;         % a number of your choice
    Lc_over_x = 4;          % larger value means more resolution
    Lc = 1.377*Gprime*Dc/(b*sigma);  %m
    Lamb = (Gprime*Dc)./((b-a).*sigma);
    Ln = (1/pi)*(Gprime*Dc)./(b.*sigma).*(b./(b-a)).^2; %nucleation length (Rubin & Ampuero 2005)
    % Yudong: 2026-06-04
    % min_fault_length = Ln*L_over_Ln; %m
    min_fault_length = 1; %m
    dx = min(Lc, Lamb) / Lc_over_x; % resolution (m)

    % coordinates
    o.topo = wl.topo;
    o.vtx = wl.vtx;

    % calculate the center coordinates 
    o.nelem = size(o.topo, 1);
    o.xc = (o.vtx(o.topo(:,1)+1,1) + o.vtx(o.topo(:,2)+1,1) + o.vtx(o.topo(:,3)+1,1))/3;
    o.yc = (o.vtx(o.topo(:,1)+1,2) + o.vtx(o.topo(:,2)+1,2) + o.vtx(o.topo(:,3)+1,2))/3;
    o.zc = (o.vtx(o.topo(:,1)+1,3) + o.vtx(o.topo(:,2)+1,3) + o.vtx(o.topo(:,3)+1,3))/3;
    o.xc = o.xc';
    o.yc = o.yc';
    o.zc = o.zc';
    one = ones(size(o.xc));

    % filling in the parameters for elements
    o.mu0 = o.mu_0*one;
    o.d_c = Dc*one;
    o.a = a*one;
    o.b = b*one;
    o.s_normal = sigma*one;
    
    o.kes=0;  %minimum normal stress (for rough fault)
   
   % Green's functions         
   % temporarily use o.mu to fix. TODO: change the matrix
    s.hm.cc{1,1} = o.mu*wl.Gs;
    if (o.use_Gn) s.hm.cc{1,2} = -o.mu*wl.Gn; end   % use negative 

    % calculate H-matrices from Greens functions matrices
    s.hm.cc{1,1} = gf2hm(s.hm.cc{1,1}, [o.hm_filename '_comp11.hm'], ...
        o.xc, o.yc, o.zc, run_cmd, o.hmmvp_mode);
    if (o.use_Gn) 
            s.hm.cc{1,2} = gf2hm(s.hm.cc{1,2}, [o.hm_filename '_comp1n.hm'], ...
        o.xc, o.yc, o.zc, run_cmd, o.hmmvp_mode);
    end
  
    % loading
    o.load = 'resolve'; % 'bc', 'uniload', 'uniload2', and 'resolve' 
    % taudot = 0; % shear stress rate [Pa/s]
    % sigmadot = 0; % normal stress rate [Pa/s]
    % Sij.sigmaD = 0.0160;   % resolved field of stress rate (e.g. rough fault with tectonic loading)
    % Sij.Psi = -45 + strike_angle;

    % Yudong: 2025-07-14, 3D loading 
    Sij.matrix6 = [0, -0.008, 0, 0, 0, 0]';
    
    switch o.load
        case 'uniload'
            o.uniload = 1;
            o.taudot = taudot;
            % Yudong May 15 2025
            if o.use_Gn
                o.IncludeNormal = true;
                o.sigmadot = 0; % consider Gn but not sigmadot (normal stressing loading)
            end
        case 'uniload2' %FIXME implement this properly in fdra2c
            o.uniload = 1;
            o.IncludeNormal = true;
            o.taudot = taudot;
            % Yudong May 6 2025 testing  
            o.sigmadot = -sigmadot;
            % o.sigmadot = 0;
        % Yudong May 29 2025
        % add tectonic loading (shear and normal stress, may be different for each fault patch)
        case 'resolve'
            o.uniload = 1;
            o.IncludeNormal = true;
            % add variable loading rate for shear and normal stress (taudot and sigmadot)
            cstress.Sij = Sij;
            cstress.topo=o.topo;
            cstress.vtx=o.vtx;

            [tau,sig]=resolve_3D_tri(cstress);
               
            o.taudot = tau';
            o.sigmadot = -sig'; % use negative 
        case 'bc'
            % Yudong Apr 16 2025
            o.uniload = 0;
            % Yudong May 15 2025
            if o.use_Gn
                o.IncludeNormal = true;
            end
    end %switch load

    if o.uniload
        o.v_creep = 1.0;
        if length(o.taudot)==1 
           o.tdot_o_vcreep=(o.taudot/o.v_creep)*one;
        else 
           o.tdot_o_vcreep=(o.taudot(:)/o.v_creep);
        end
        % Yudong May 5 2025
        % add normal stress loading rate
        if o.IncludeNormal 
           if length(o.sigmadot)==1 
             o.ndot_o_vcreep=(o.sigmadot/o.v_creep)*one;
           else 
             o.ndot_o_vcreep=(o.sigmadot(:)/o.v_creep);
           end
        end
   
    else % case 'bc'
        o.v_creep = Vpl_creep;
        % Yudong Apr 16 2025, read in bc
        % Green's function relating stressing rates to BC slip rates
        o.hm_bc = Gs(:,end); % Gload
       
        % Yudong May 14 2025
        if o.use_Gn
           o.hm_bc_n = -Gn(:,end); % Gload_n, use negative
        end
   
    end

    % start and end time
    o.ti = 0;
    % o.tf = 8e8;
    o.tf = 2e9;

    % initial conditions 
    o.slip_init = one*0;
    o.v_init = max(Vpl_creep, 1e-11)*one;
    o.chi_init = one*0;
    if o.ncomp==2
        o.v_init = [o.v_init; o.vzero*one];
    end

    % -------- run commands, write output files ----------- %

    % hard link H-matrices with .hmat extension (compatible with fdra2c)
    % compXY indicates source and receiver components (1=along strike, 2=along dip, n=normal).
    % we're assuming no opening so there is no normal source.
    if (run_cmd)
        % create symbolic links to H-matrices
        system(['ln -s ' o.hm_filename '_comp11.hm ' o.hm_filename '_comp11.hmat']);
        if o.use_Gn system(['ln -s ' o.hm_filename '_comp1n.hm ' o.hm_filename '_comp1n.hmat']); end
        if o.ncomp==2
            system(['ln -s ' o.hm_filename '_comp12.hm ' o.hm_filename '_comp12.hmat']);
            system(['ln -s ' o.hm_filename '_comp21.hm ' o.hm_filename '_comp21.hmat']);
            system(['ln -s ' o.hm_filename '_comp22.hm ' o.hm_filename '_comp22.hmat']);
            if o.use_Gn system(['ln -s ' o.hm_filename '_comp1n.hm ' o.hm_filename '_comp2n.hmat']); end
        end
    end

    % write out readme file
    e.description=['3D_triangular_fault'];
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


end % function simul_setup



function addpaths ()
    e = env();
    addpath([e.work_dir '/fdra2c/matlab']);
end

function o = set_defaults();

    e=env();
    o.ncomp = 1; %number of slip components. 1 for 1D, 1 or 2 for 2D.

    % Yudong: 2025-06-25 have green's functions for normal stresses now in fdra2c!
    o.use_Gn = 0; %default: do not include greens functions for normal stresses
    o.dec=[e.dec 'no_Gn']; %add suffix to filenames to indicate whether normal stresses are included.
    o.prefix=e.prefix;
    o.dir = e.scratch_dir;
    
    o.hmmvp_mode = 'omp'; %parallelization mode for hmmvp.
    o.stress_fn = 'h-matrix';
    
    o.tol = 1e-9; % tolerance
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

end % function set_defaults

function [tau, sig] = resolve_sigmaD(c)
    % Functions to resolve stresses (e.g. Schultz, Geological Fracture Mechanics, eq.3.4)
    % Assume mean stress = 0.
    % sigmaD is the angle the fault makes to sigma1. 
    % Sign conventions: CCW=+ve; left-lateral=+ve; compressive=+ve.
    
    if c.multiseg
       dx=[];  dy=[];
       for n=1:c.nseg
          dx = [dx diff(c.xs{n})];
          dy = [dy diff(c.ys{n})];
       end
    else
       dx = diff(c.xs);
       dy = diff(c.ys);
    end
    
    theta= atan2d(dy, dx); %angle from fault plane
    tau= -c.Sij.sigmaD/2*sind(2*(c.Sij.Psi+theta)).*c.sense;
    sig= -c.Sij.sigmaD/2*cosd(2*(c.Sij.Psi+theta));
    
    tau=tau(:);
    sig=sig(:);
end % function resolve_sigmaD

function [tau, sig] = resolve_3D(c)
    % Yudong: 2025-07-14, 3D loading 
    tau = zeros(size(c.elr,2),1);
    sig = zeros(size(c.elr,2),1);
    for i = 1:size(c.elr,2)
        along=c.elr(6,i)*[sind(c.elr(5,i)),cosd(c.elr(5,i)),0]';
        normal=[cosd(-c.elr(5,i))*sind(c.elr(4,i)),sind(-c.elr(5,i))*sind(c.elr(4,i)),cosd(c.elr(4,i))]';
        [ss ns] = Stresses(c.Sij.matrix6,along,normal);
        tau(i) = ss;
        sig(i) = ns;
    end
end % function resolve_3D

% Yudong:2026-06-12, calculate the resolved stressing for triangular elements.
% check the geometry
function [tau, sig] = resolve_3D_tri(c) 
    tau = zeros(size(c.topo,1), 1);  
    sig = zeros(size(c.topo,1), 1);
    
    for i = 1:size(c.topo,1)  
        % get the coordinates of the three vertices of the triangle
        v1 = c.vtx(c.topo(i,1)+1, :);  
        v2 = c.vtx(c.topo(i,2)+1, :);
        v3 = c.vtx(c.topo(i,3)+1, :);
        
        % calculate the edge vectors of the triangle
        edge1 = v2 - v1;  
        edge2 = v3 - v1;  
        
        % calculate the normal vector of the triangle
        % normal = cross(edge1, edge2);
        normal = cross(edge2, edge1);  
        normal = normal / norm(normal);  
        
        
        vert = [0, 0, 1];  % veritcal vector
        
        % along strike direction.
        % along = cross(normal,vert);
        along = cross(vert,normal);
        along = along / norm(along); 
        
        % stress
        % [ss ns] = Stresses(c.Sij.matrix6, along', -normal');
        [ss ns] = Stresses(c.Sij.matrix6, along', normal');
        tau(i) = ss;
        sig(i) = ns;
    end
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
       %Yudong May 6 2025 testing
       % c.tol = 1e-6;
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
    
end  % function gf2hm

% -------- Write key-value file for fdra2c ----------- %
function cf = simul_write_fdra_kvf (o)
    e=env();
    cf=o;
    cf.kvf = [o.save_filename '_f'];
    cf.cmd = sprintf(['mpirun %s/fdra2c/bin/fdra %s.kvf;\n'], e.work_dir, cf.kvf);
    kvf('Write',[cf.kvf '.kvf'], o, true);
end % function simul_write_fdra_kvf

function [Gs Gn] = GF_disloc3d_ss(els,elr,mu,nu)
    % compute the Green's functions for disloc3d for strike slip case
    % els: sources (x, y, z, length, width, dip, strike,sense)
    % elr: receivers (x, y, z, dip, strike,sense), it means (east, north, up), z is negative 
    % mu: shear modulus [Pa]
    % nu: Poisson's ratio
    % Gs: shear Green's functions
    % Gn: normal Green's functions
    % left_lateral: 1 for left-lateral, -1 for right-lateral
    % Yudong: 2025-06-30
    Gs = zeros(size(elr,2), size(els,2));
    Gn = zeros(size(elr,2), size(els,2));
    obs1 = elr(1:3,:);
    for j = 1:size(els,2)
        % model parameters: mdl = [length width depth dip strike east north ss ds op]';
        mdl_para=[els(4,j), els(5,j), -els(3,j), els(6,j), els(7,j), els(1,j), els(2,j), els(8,j), 0, 0]';
        [Uf Df Sf flag] = disloc3d(mdl_para,obs1,mu,nu);
        for i = 1:size(elr,2)
            along=elr(6,i)*[sind(elr(5,i)),cosd(elr(5,i)),0]';
            normal=[cosd(-elr(5,i))*sind(elr(4,i)),sind(-elr(5,i))*sind(elr(4,i)),cosd(elr(4,i))]';
            [ss ns] = Stresses(Sf(:,i),along,normal);
            Gs(i,j) = ss;
            Gn(i,j) = ns;
        end
    end
end

% Yudong: 2025-07-14, vertical integral
function [Gs Gn] = GF_disloc3d_ss_vertical_integral(els,elr,mu,nu)
    % compute the Green's functions for disloc3d for strike slip case
    % els: sources (x, y, z, length, width, dip, strike,sense)
    % elr: receivers (x, y, z, dip, strike,sense), it means (east, north, up), z is negative 
    % mu: shear modulus [Pa]
    % nu: Poisson's ratio
    % Gs: shear Green's functions
    % Gn: normal Green's functions
    % left_lateral: 1 for left-lateral, -1 for right-lateral
    % Yudong: 2025-06-30
    Gs = zeros(size(elr,2), size(els,2));
    Gn = zeros(size(elr,2), size(els,2));
    obs1 = elr(1:3,:);
    for j = 1:size(els,2)
        % model parameters: mdl = [length width depth dip strike east north ss ds op]';
        % mdl_para=[els(4,j), els(5,j), -els(3,j), els(6,j), els(7,j), els(1,j), els(2,j), els(8,j), 0, 0]';
        % Yudong comment: 2025-07-14, nk need to be an odd number, don't know why
        nk = 2*ceil(els(5,j)/els(4,j)/2)+1;
        % nk = 9; % for testing
        Sf_vertical_integral = zeros(6,size(elr,2));
        for k = 1:nk
            wid_tmp = els(5,j)/nk;
            dep_tmp = -els(3,j) - (-1+k)*wid_tmp*sind(els(6,j));
            est_tmp = els(1,j) - (-1+k)*wid_tmp*cosd(els(6,j))*cosd(els(7,j));
            nth_tmp = els(2,j) + (-1+k)*wid_tmp*cosd(els(6,j))*sind(els(7,j));
            % normalized slip
            s_tmp = sqrt(1-((-1/2+k)*wid_tmp-els(5,j)/2)^2/(els(5,j)/2)^2); 
            mdl_para=[els(4,j), wid_tmp, dep_tmp, els(6,j), els(7,j), est_tmp, nth_tmp, s_tmp*els(8,j), 0, 0]';
            [Uf Df Sf flag] = disloc3d(mdl_para,obs1,mu,nu);
            Sf_vertical_integral = Sf_vertical_integral + Sf;
        end
        for i = 1:size(elr,2)
            along=elr(6,i)*[sind(elr(5,i)),cosd(elr(5,i)),0]';
            normal=[cosd(-elr(5,i))*sind(elr(4,i)),sind(-elr(5,i))*sind(elr(4,i)),cosd(elr(4,i))]';
            [ss ns] = Stresses(Sf_vertical_integral(:,i),along,normal);
            Gs(i,j) = ss;
            Gn(i,j) = ns;
        end
    end
end

function [ss ns] = Stresses(S,along,normal)
    % From S, the output of disloc3d, compute the shear and normal stresses. along
    % and normal are vectors pointing along and normal to the fault.
    n = size(S,2);
    ss = zeros(n,1);
    ns = zeros(n,1);
    % Vectorized code that is equivalent to
    %   sigma = [Sxx(i) Sxy(i) Sxz(i)
    %            Sxy(i) Syy(i) Syz(i)
    %            Sxz(i) Syz(i) Szz(i)];
    %   ss(i) = along' *sigma*normal;
    %   ns(i) = normal'*sigma*normal;
    normal = normal(:);
    S = S([1 2 3 2 4 5 3 5 6],:);
    nml = repmat(normal(repmat(1:3,1,3)),1,n);
    a = S.*nml;
    a = [sum(a(1:3,:))
        sum(a(4:6,:))
        sum(a(7:9,:))];
    ss = along(:)'*a;
    ns = normal(:)'*a;
end % function Stresses

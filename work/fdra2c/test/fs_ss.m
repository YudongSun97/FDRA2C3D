% Input script to fdra(). Copy this file and modify the values to set up
% different simulations. Call as fdra('fs').

% Evolution law: aging (1), logarithmic (2), Perrin (3)
evolution = 2;

% Dilatancy law
dilatancy_law = 0;
dlt_phi0 = 0.05;

% Constants
year2sec = 365.25*24*60*60;
sec2year = 1/year2sec; 

% Constitutive parameters
mu_0    = 0.60;              % nominal friction coefficient
epsilon = 5.0e-5;
beta    = 6.0e-11;           % compressibility (1/Pa) = phi*(beta_f + beta_phi)
Chyd    = 1.0d-6;            % Hydraulic Diffusivity [m^2/s]
Cth     = 1.0e-6;            % Thermal Diffusivity [m^2/s]
G       = 3.0d10;            % Rigidity [Pa]
poisson = 0.25;              % Poisson ratio 
Gprime  = G/(1-poisson);     % G/(1-nu) for plane strain; G for anti-plane
V0      = 1.d-6;             % Reference velocity [m/s]
Vs      = 3.7d3;             % S wave velocity [m/sec]
eta     = 0.5*Gprime/Vs;     % radiation damping parameter
h_c     = 1e-3;
h       = 1e-3;              % shear zone thickness [m]
hfac = 1;
h = hfac*h;
% Temperature
rho     = 2600;              % Total density [kg/m^3]
c_heat  = 1100;              % Heat Capacity [N-m/kg/deg C]
Lambda  = 8e5;  % 0 to turn off. Thermal Pressurization Coeff  [Pa /deg C]
% Linker-Dieterich
LD_alpha = 0.3;
% Boundary conditions
Vpl_creep = 0.04*sec2year;   % secular slip in a deep region [cm/yr]    

% spring-slider model
nrep = 40;
o = ones(nrep,1);
elasticity = 'ss';
mu_0 = mu_0*o;
b = 0.018*o;
a = 0.9*b;
Dc = 1e-4*o;
s_normal = 100e6*o;
p_inf = 0;
k_crit = (s_normal - p_inf).*(b - a)./Dc;
ss_k = 0.01*k_crit;

x = 0*o;
Nflt = nrep;
N_cell = Nflt;

Prem   = p_inf;
one    = ones(Nflt,1);
rho    = rho*one;
c_heat = c_heat*one;

% Mesh parameters normal to fault
Ny = 150; % Number of spatial cells
z_c = 1e-5;
% Chyd and Cth can vary along the fault. Set the value for each x point.
Chyd = Chyd*ones(Nflt,1);
Cth  = Cth*ones(Nflt,1);
y_inf = 10*sqrt(Chyd.*30*year2sec); % Diffusion length for 30 year cycle time

% Finite-width shear zone
mdsbi_mode = 1;
sz_fw = 0;
if(sz_fw)
  % Number of cells to dedicate to the shear zone
  sz_Nh = ceil(0.5*Ny);
  Ny = [sz_Nh Ny];
  sz_Chyd = min(Chyd)*one;
  sz_Cth  = min(Cth)*one;
  % Falloff parameter
  sz_fo = sqrt(2)/min(h);
  %sz_fo = 1e2;
end
sNy = sum(Ny);

% Initial conditions
v_init = 10*Vpl_creep*ones(size(a));
psi_init = log(v_init/V0);
% Initial State    
chi_init = zeros(size(psi_init));
foo = (p_inf'*ones(1,sNy))';
p_init  = foo(:);
T_init = zeros(size(p_init));

% Time interval
ts   = 0;
tend = 1e9;

% Profile locations
profs = 1;

% Flash heating
flash_heating = 0;
flash_fw      = 0.13;
flash_tauc    = 3e9;
flash_Tw      = 300;
flash_D       = 5e-6; % Asperity diameter

v1 = 1e2*o;
v2 = 1e-6*o;
use_vcutoff = 1;

rsUse = 0;

Chyd = inf;
Lambda = 0;
LD_alpha = 0;
no_dilatancy = 1;
odepair = '23';

% Solver parameters
% Uncomment one or both of these to turn off:
use_Gn    = 1;          % Variable normal stress
kludge    = 1;          % Protect against p > sigma in shear zone
use_mixed = 1;          % mixed DAE/ODE time stepping
pure_ode = 1; fprintf(1,'----> PURE ODE\n');
use_compressed_bem = 1; % compressed BEM (small err, large speed up)
use_nthreads = 4;       % Number of OpenMP threads to use
monitor_threads = 0;

Nsavepoints = 10;        % How often to save data in solver time steps
conType = 'loo';      % Convergence tolerance: loose, moderate, tight
dispFig = 0;            % Plot results as solver works
dispText = 1000;

% Solver parameters
t = ''; l = ''; vs = ''; mixed = ''; ngn = ''; k = ''; fwsz = ''; fh = '';
rsStr = '';
if(LD_alpha)      l     = 'l'; end
if(Lambda)        t     = 't'; end
if(use_mixed)     mixed = 'm'; end
if(~use_Gn)       ngn   = 'nGn'; end
if(sz_fw)         fwsz  = 'fw'; end
if(kludge)        k     = 'k'; end
if(flash_heating) fh    = 'fh'; end
saveFn = '/scratch/ambrad/fdra2c/mcode';

% Protect against overwriting existing simulation data
if(~isempty(dir([saveFn '.mat'])))
%  error('File exists, my friend!');
end

print_debug = 0;

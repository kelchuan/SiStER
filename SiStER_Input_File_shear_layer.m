% SiStER_Input_File


% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=100; % max number of time iterations
dt_out=10; % output files every "dt_out" iterations


% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=15;
ysize=3;
% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y
GRID.dx(1)=0.1;
GRID.x(1)=15;
GRID.dy(1)=0.1;
GRID.y(1)=3;



% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=6; % number of markers in the smallest quadrant
Mquad_crit=3; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphase=1; % number of phases

% phase 1 
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=3;



% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1
MAT(1).phase=1;
% density parameters
MAT(1).rho0=3000;
MAT(1).alpha=0;
% elasticity 
MAT(1).G=1e10;
% diffusion creep parameters
MAT(1).pre_diff=.5/1e20;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=.5/1e20;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).Cmax=40e6;
MAT(1).Cmin=0.01e6;
MAT(1).ecrit=0.1;





% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=0; % elasticity on (1) or off (0)
PARAMS.YNPlas=0; % plasticity on (1) or off (0)
PARAMS.tau_heal=1e11; % healing time for plasticity (s)
PARAMS.gx=0; % gravity along x
PARAMS.gy=0; % gravity along y
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=1e25; % maximum viscosity
PARAMS.etamin=1e18; % minimum viscosity
PARAMS.Tsolve=0; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.a0=0;
PARAMS.a1=0;
PARAMS.a2=0;
PARAMS.a3=1000/(30e3)^3;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=0;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(1).rho0; 
PARAMS.kref=3;
PARAMS.cpref=1000;


% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-8; % diffusivity of topography (m^2/s)

% Picard iterations
PARAMS.Npicard_min=10; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=100; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-9;
PARAMS.pitswitch=0; % number of Picard iterations at which the solver switches to quasi-Newton


% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=1e9; % pressure in the top-left corner of the domain (anchor point)

%velocity function, consistent with strain rate for shear
y=0;
for i=1:length(GRID.y)
  n=round((GRID.y(i)-y(end))./GRID.dy(i));
  dd=(GRID.y(i)-y(end))./n;
  y=[y, y(end)+(1:n).*dd];
end
n=round((ysize-y(end))./GRID.dy(end));
dd=(ysize-y(end))./n;
y=[y, y(end)+(1:n).*dd];

PARAMS.boxheight = y;
PARAMS.strain_rate_shear = 1e-15;
PARAMS.vmax_shear = PARAMS.strain_rate_shear*ysize;   
PARAMS.velocity_profile=@(position) (PARAMS.vmax_shear*(position/ysize)); 
% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant)
% 3/ value of normal velocity 

BC.top=[0 0 0];
BC.bot=[0 0 0 0];
BC.left=[1 0 0];
BC.right=[1 0 0];
BC.bot_profile=-(.01/365.25/24/3600)*ones(1,481).*(1-exp(-([1:481]-1)/10));
BC.vshear = PARAMS.vmax_shear;
PARAMS.BalanceStickyLayer=0; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 2=Neumann)
% 2/ value
BCtherm.top=[1 550];
BCtherm.bot=[1 550];
BCtherm.left=[2 0];
BCtherm.right=[2 0];

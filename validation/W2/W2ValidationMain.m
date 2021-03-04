%---------------------------- Header -------------------------------------%
% Title: W2 validation input file
% Author: J. Davis
% Date: 12-15-20
% Last Update: 12-15-20
% Description: Validation input file (Main.m) for W2 experimental and
% numerical WAMIT/WEC-Sim data.
%
% Sources:
%   (1) Cite M. Choiniere here
%-------------------------------------------------------------------------%

% Validation data:
Val.Source1 = 'W2 WAMIT';    load('.//W2_hydro.mat');    Val.Data1 = W2_hydro;
% % Val.Source2 = 'W2 WEC-Sim';  load('.//W2_MCR_free.mat'); Val.Data2 = W2_MCR;
Val.Source2 = 'W2 WEC-Sim';  load('.//W2_MCR_0degreg_interp2.mat'); Val.Data2 = W2_MCR_0degreg ;
Val.Source3 = 'W2 0deg exp'; load('.//W2_0deg_exp.mat'); Val.Data3 = W2_0deg_exp;

% Constants;
env.g       = 9.81;             % Acceleration of gravity [m/s^2]
env.rho     = 1025;             % Density of sea water [kg/m^3]

% Environment parameters
env.omega   = Val.Data1.w;      % Wave angular frequencies [rad/s]
env.T       = 2*pi./env.omega;  % Wave periods [s]
env.Aw      = .06;              % Wave amplitudes [m]
env.defwavesteepness = 0;       % Parameterize wave amplitudes based on steepness?
    env.steepness = 0.0049;     % Wave steepness, H/L [-]
env.h       = 4.5;              % Water depth [m]

% Use T and A from experiments:
% env.omega   = Val.Data2.w;        % Wave angular frequencies [rad/s]
% env.T       = 2*pi./env.omega;    % Wave periods [s]
% env.Aw      = Val.Data2.H/2;      % Wave amplitudes [m]
% env.h       = 4.5;                % Water depth [m]


% OSWEC dimensions:
w = 0.94;      % Width [m]
c = 3.848;     % Dist from bottom [m]
t = 0.095;     % Thickness [m]
ht = env.h-c;  % Height [m]

% Assign to param object (only neccesary for single dimension runs):
param.w_list = w;   % Width [m]
param.c_list = c;   % Distance from bottom [m]; c = water depth - device height
param.w2tr   = w/t;  % Width to thickness ratio [m/m]

% Mass properties
body.defbodyprop  = true;                           % Use user-defined body properties?
body.prop.rho_m   = env.rho/2;                      % Structural density [kg/m^3]
body.prop.mass    = 25.7;                           % Device mass [kg]
body.prop.I55     = body.prop.rho_m/3 ...           % Pitch mass moment of inertia [kg-m^2]
    *w*t*ht^3*(1+(t/(2*ht))^2);
  
body.prop.rb      = ht-0.294;                       % Buoyant torque arm [m]
body.prop.rg      = ht-0.364;                       % Gravity torque arm [m]
body.hydro.Badd   = 7.5;                           % Additional damping [kg-m^2/s]
body.hydro.C55    = (env.rho*w*t*ht*body.prop.rb...
            - body.prop.mass*body.prop.rg)*env.g;   % Device's buoyant torque calculated from buoyancy [N-m]
body.paramBodyProp= true;
% body.hydro.C55     = 27; % This is the number from wecSim, about CoG

% PTO properties:
body.pto.ctrltype = 'free';             % PTO control scheme: 'free','damper','CCC'
body.pto.PTO_eff = 1.0;                 % PTO efficiency [-]
    % CCC control settings:
    body.pto.pitch_constraint = true;   % Apply pitch constraints? (Only for CCC)
    body.pto.xi5_max = 30*pi/180;       % Max pitch amplitude [rad]
    % Damper control settings:
    body.pto.Cg_constant = 0;           % Constant PTO restoring coefff [Kg-m^2/s^2]
    solver.calculateACE = false;        % Calculate ACE? (Only works for damper)
    % Free settings:
    solver.calculateHydroEff = false;   % Calculate hydrodynamic efficiency?
    
% Foundation properties:
fdn.prop.rho = 2450;        % Density of reinforced concrete [kg/m^3]
fdn.prop.Vratio = 1;        % Ratio of foundation volume to OSWEC volume
fdn.paramFdnProp = false;

% Convergent factors for solving dispersion (wn) and numerical scheme:
solver.n = 4;               % Dispersion equation - number of frequencies to keep
solver.nmax = 4;            % Collocation scheme for Chebyshev problem or Mathieu solutions

% Parallel processing:
solver.parallel = 1;        % 1- On; Otherwise - Off. For Renzi-Diaz' model only.

% Post-processing file:
% PostProcessFile = './/validation//W2//W2PostProcess.m';
PostProcessFile = './/validation//W2//W2PostProcessEWTEC21.m';
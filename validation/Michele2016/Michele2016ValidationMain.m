%---------------------------- Header -------------------------------------%
% Title: W2 validation input file
% Author: J. Davis
% Date: 12-15-20
% Last Update: 12-15-20
% Description: Validation input file (Main.m) for W2 experimental and
% numerical WAMIT/WEC-Sim data.
%
% Sources:
%   (1) N. Tom, Y.-H. Yu, and A. Wright, "Balancing Power Absorption and
%       Structural Loading for a Novel Fixed-Bottom Wave Energy Converter
%       with Nonideal Power Take-Off in Regular Waves: Preprint," p. 6, 2017.
%-------------------------------------------------------------------------%

% Constants;
env.g       = 9.81;             % Acceleration of gravity [m/s^2]
env.rho     = 1025;             % Density of sea water [kg/m^3]

% Environment parameters
env.omega   = 0.25:0.01:4.0;    % Wave angular frequencies [rad/s]
env.T       = 2*pi./env.omega;  % Wave periods [s]
env.Aw      = 1;                % Wave amplitudes [m]
env.h       = 5;                % Water depth [m]

% OSWEC dimensions:
w = 3;         % Width [m]               
c = 0;         % Dist from bottom [m]            
t = 0.25;      % Thickness [m]          
ht = env.h-c;  % Height [m]
%
param.w_list = w;   % Width [m]
param.c_list = 0;   % Distance from bottom [m]; c = water depth - device height 
param.w2tr   = w/t;  % Width to thickness ratio [m/m]

% Mass properties
body.defbodyprop = true;
body.prop.rho_m   = 36*10^3/72;         % Structural density [kg/m^3]
body.prop.mass    = w*t*ht*rho_m;       % Device mass [kg]
body.prop.rb      = ht/2;               % Buoyant torque arm [m]
body.prop.rg      = ht/2;               % Gravity torque arm [m]
body.prop.I55     = 5E4;                % Pitch mass moment of inertia [kg-m^2]
body.hydro.C55    = 7E5;                % Device's buoyant torque calculated from buoyancy [N-m]
        
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

% Convergent factors for solving dispersion (wn) and numerical scheme;
% Generally, convergence can be achieved with n = 4 and nmax = 4; 
solver.n = 4;               % Dispersion equation - number of frequencies to keep
solver.nmax = 4;            % Collocation scheme for Chebyshev problem or Mathieu solutions

% Parallel processing;
solver.parallel = 1;        % 1- On; Otherwise - Off. For Renzi-Diaz' model only.

% HAVENT PREPARED THESE YET:
% Validation data:
Val.Source1  = '';    load('');    Val.Data1 = ;

% Post-processing file:
PostProcessFile = './/validation//Michele2016//';


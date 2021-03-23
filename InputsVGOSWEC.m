%---------------------------- Header -------------------------------------%
% Title: InputsVGOSWEC.m
% Date: 3-10-21
% Last Update: 3-10-21
% Description: Input file to Main.m for VGOSWEC scale model sizing
%-------------------------------------------------------------------------%

% Constants 
env.g = 9.81;               % Acceleration of gravity [m/s^2]
env.rho = 1025;             % Density of sea water [kg/m^3]

% Environment parameters
env.h = 1;                   % Water depth [m]
env.T = 0.5:0.01:3;          % Wave periods [s]
env.omega = 2*pi./env.T;     % Wave angular frequencies [rad/s]
% env.omega = 0.1:0.5:15;     % Wave angular frequencies [rad/s]
% env.T = 2*pi./env.omega;    % Wave periods [s]

env.defspectrum = 0;
    env.Tp = 9.86;
    env.Hs = 2.64;
env.Aw = 0.1;               % Wave amplitudes [m]
env.defwavesteepness = 1;   % Parameterize wave amplitudes based on steepness?
    env.steepness = 0.0049;     % Wave steepness, H/L [-]
    env.steepness = 0.01;     % Wave steepness, H/L [-]
    
% OSWEC parameterizations:
tank_width = 1.2; %[m]
% param.w_list = 1/4*tank_width:.1:1/2*tank_width;% Device width [m]
param.w_list = 0.4; % [m]
%param.c_list = 0:.1:1/2*env.h; % Distance from bottom [m]; c = water depth - device height
param.c_list = 0.5; % Distance from bottom [m]; c = water depth - device height
    

% param.w_list = 0.3:0.1:0.8; % [m]
% param.c_list = 1/2*env.h; % Distance from bottom [m]; c = water depth - device height

% param.rho_m_list = [300 400 500 600 700 800 900 1000];
% param.w_list = 0.4; % [m]
% param.c_list = 0.5*ones(size(param.rho_m_list)); % Distance from bottom [m]; c = water depth - device height

% param.w2tr = [8  7  6  5  4].';
% param.w_list = 0.4; % [m]
% param.c_list = 0.5*ones(size(param.w2tr)); % Distance from bottom [m]; c = water depth - device height


param.w2tr = 8;         % Width to thickness ratio [m/m]
body.prop.rho_m = env.rho/2; % Structural density [kg/m^3]
body.paramBodyProp= true;
body.Bmech = 2; % Mechanical damping in the system [kg-m^2/s] 

% PTO properties:
body.pto.ctrltype = 'free';           % PTO control scheme: 'free','damper','CCC'
body.pto.PTO_eff = 1.0;                 % PTO efficiency [-]
    % CCC control settings:
    body.pto.pitch_constraint = false;  % Apply pitch constraints? (Only for CCC)
    body.pto.xi5_max = 30*pi/180;       % Max pitch amplitude [rad]
    % Damper control settings:
    body.pto.Cg_constant = 0;           % Constant PTO restoring coefff [Kg-m^2/s^2]
    solver.calculateACE = false;        % Calculate ACE? (Only works for unconstrained)
    % Free settings:
    solver.calculateHydroEff = false;
    
% Foundation properties:
fdn.paramFdnProp = true;
    fdn.prop.rho = 2450;        % Density of reinforced concrete [kg/m^3] 
    fdn.prop.SigYield = 70E6;   
    % fdn.prop.Vratio = 1;        % Ratio of foundation volume to OSWEC volume
    fdn.prop.SF = 1.75;
    fdn.prop.gamma = 0.75; % Ratio of annular cross section outer to inner diameter
fdn.prop.mass = [];
%%% IN THE FUTURE, INCLUDE MATERIAL SELECTION OPTION

% Convergent factors for solving dispersion (wn) and numerical scheme;
solver.n = 4;               % Dispersion equation - number of frequencies to keep
solver.nmax = 4;            % Collocation scheme for Chebyshev problem or Mathieu solutions

% Parallel processing;
solver.parallel = 1;        % 1- On; Otherwise - Off. For Renzi-Diaz' model only.

% Post-processing file:
PostProcessFile = 'PostProcessVGOSWEC.m';
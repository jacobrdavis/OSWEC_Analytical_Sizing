%---------------------------- Header -------------------------------------%
% Title: Generic OSWEC Analytical Sizing input file
% Author: J. Davis
% Date: 7-23-20
% Last Update: 12-15-20
% Description: This is an analytical approach to optimizing the dimensions
%   of a bottom-hinged rectangular paddle wave energy converter based on
%   time-average power. This script utilizes the "Semi-analytical and 
%   analytical method to evaluate hydrodynamic coefficients (added mass 
%   and damping) of an surface piercing OSWEC" tool developed by Jessica 
%   Nguyen, PhD. 
%
% Sources:
%   (1) N. Nguyen, "Semi-analytical and analytical method to evaluate
%       hydrodynamic coefficients (added mass and damping) of an surface
%       piercing OWEC," MATLAB Script, July 2020
%
%   (2) S. Michele, P. Sammarco, and M. d'Errico, "Theory of the
%       synchronous motion of an array of floating flap gates oscillating 
%       wave surge converter," Proc. R. Soc. A, vol. 472, no. 2192, p. 
%       20160174, Aug. 2016, doi: 10.1098/rspa.2016.0174.
%
%   (3) N. Tom, Y.-H. Yu, and A. Wright, "Balancing Power Absorption and
%       Structural Loading for a Novel Fixed-Bottom Wave Energy Converter
%       with Nonideal Power Take-Off in Regular Waves: Preprint," p. 6, 2017.
%
%   (4) J. N. Newman, "The Exciting Forces on Fixed Bodies in Waves,"
%       Journal of Ship Research, The Society of Naval Architects and Marine
%       Engineers, vol. 6, no. 4, Dec. 1962.
%-------------------------------------------------------------------------%
close all
clear env param body fdn solver out
%% INPUTS AND SETTINGS

% Constants 
env.g = 9.81;               % Acceleration of gravity [m/s^2]
env.rho = 1025;             % Density of sea water [kg/m^3]

% Environment parameters
env.T = 0.5:0.01:4;          % Wave periods [s]
env.omega = 2*pi./env.T;     % Wave angular frequencies [rad/s]
% env.omega = 0.1:0.5:15;     % Wave angular frequencies [rad/s]
% env.T = 2*pi./env.omega;    % Wave periods [s]
env.Aw = 0.06;               % Wave amplitudes [m]
env.Asc = 1;                 % Wave amplitude scale [-]
env.h = 1;                   % Water depth [m]

% OSWEC parameterizations:
param.w_list = 0.5;          % Device width [m]
param.c_list = 0:0.1:0.5;    % Distance from bottom [m]; c = water depth - device height
param.w2tr = 20/1.5;         % Width to thickness ratio [m/m]
body.prop.rho_m = env.rho/2; % Structural density [kg/m^3]
body.defbodyprop = false;

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
    solver.calculateHydroEff = true;
    
% Foundation properties:
fdn.prop.rho = 2450;        % Density of reinforced concrete [kg/m^3] 
fdn.prop.Vratio = 1;        % Ratio of foundation volume to OSWEC volume

% Convergent factors for solving dispersion (wn) and numerical scheme;
solver.n = 4;               % Dispersion equation - number of frequencies to keep
solver.nmax = 4;            % Collocation scheme for Chebyshev problem or Mathieu solutions

% Parallel processing;
solver.parallel = 1;        % 1- On; Otherwise - Off. For Renzi-Diaz' model only.

Post-processing file:
PostProcessFile = 'PostProcess.m';

%% Test Case: W2 Experiments
% run('.//validation//W2//W2ValidationMain.m')

%% Test Case: Tom 2017
% run('.//validation//Tom2017//Tom2017ValidationMain.m')

%% Test Case: Michele 2016
% run('.//validation//Michele2016//Michele2016ValidationMain.m')

%% CALCULATION LOOP
[out,env,body,fdn,param,solver] = fun.OSWEC_CalcLoop(env,body,fdn,param,solver);

%% POST PROCESSING
run(PostProcessFile)

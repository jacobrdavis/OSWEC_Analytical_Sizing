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
%close all
clear env param body fdn solver out
%% INPUTS AND SETTINGS
% InputFile = 'InputsVGOSWEC.m';
InputFile = 'InputsEWTEC21.m';
% InputFile = './/validation//W2//W2ValidationMain.m';                        % Test Case: W2 Experiments
% InputFile = './/validation//Tom2017//Tom2017ValidationMain.m';              % Test Case: Tom 2017
% InputFile = './/validation//Michele2016//Michele2016ValidationMain.m';      % Test Case: Michele 2016

run(InputFile)

%% CALCULATION LOOP
[out,env,body,fdn,param,solver] = fun.OSWEC_CalcLoop(env,body,fdn,param,solver);

%% POST PROCESSING
run(PostProcessFile)
open(PostProcessFile)

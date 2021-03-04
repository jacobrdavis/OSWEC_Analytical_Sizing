function [out,env,body,fdn,param,solver] = OSWEC_CalcLoop(env,body,fdn,param,solver)
%% INITIALIZATION
nrows = length(param.c_list);
ncols = length(param.w_list);

if body.paramBodyProp == true
body.prop.I55   = zeros(nrows,ncols);
body.prop.rb    = zeros(nrows,ncols);
body.prop.rg    = zeros(nrows,ncols);
body.prop.mass  = zeros(nrows,ncols);
body.hydro.C55  = zeros(nrows,ncols);
body.hydro.Badd  = zeros(nrows,ncols);
end

fdn.prop.mass   = zeros(nrows,ncols);
fdn.dim.L   = zeros(nrows,ncols);
fdn.dim.D   = zeros(nrows,ncols);

body.hydro.mu55 = cell(nrows,ncols);
body.hydro.nu55 = cell(nrows,ncols);
body.hydro.mu15 = cell(nrows,ncols);
body.hydro.nu15 = cell(nrows,ncols);
body.hydro.F5   = cell(nrows,ncols);
body.hydro.F1   = cell(nrows,ncols);
body.hydro.X5   = cell(nrows,ncols);
body.hydro.X1   = cell(nrows,ncols);
body.pto.Cg     = cell(nrows,ncols);
body.pto.nu_g   = cell(nrows,ncols);
out.TAP         = cell(nrows,ncols);
out.CWR         = cell(nrows,ncols);
out.maxCWR      = zeros(nrows,ncols);
out.ACCW        = zeros(nrows,ncols);
out.ACETAP      = cell(nrows,ncols);
out.mtlmass     = cell(nrows,ncols);
out.CCE         = cell(nrows,ncols);
out.ACE         = cell(nrows,ncols);
out.RAO         = cell(nrows,ncols);
out.xi5         = cell(nrows,ncols);
out.Xr1         = cell(nrows,ncols);
out.Xr3         = cell(nrows,ncols);
out.Fr1         = cell(nrows,ncols);
out.Fr3         = cell(nrows,ncols);
out.maxFr1      = zeros(nrows,ncols);
out.projAreaWaveKEperW  = cell(nrows,ncols);
out.rotKE               = cell(nrows,ncols);
out.HydroEff            = cell(nrows,ncols);
out.MbFdnBaseMax   = zeros(nrows,ncols);

if size(env.omega,1)==1
    env.omega = transpose(env.omega);
    env.T = transpose(env.T);
end
if size(env.T,1)==1
    env.T = transpose(env.T);
end
if exist('env.Aw','var')
    if size(env.Aw,1)==1
        env.Aw = transpose(env.Aw);
    end
    if  length(env.Aw) == 1
        env.Aw = env.Aw*ones(size(env.omega));
    end
end
[body.dim.w, body.dim.c] = meshgrid(param.w_list,param.c_list);
%% CALCULATIONS
tic

% Parameterized device thickness and height:
body.dim.t = body.dim.w/param.w2tr;     % Thickness [m]
body.dim.ht = env.h-body.dim.c;         % Height (assumed surface piercing) [m]

% Calculate wave numbers [m^-1]
[env.k] = wave_disp(env.omega,env.h,env.g,0.001);
env.lambda = 2*pi./env.k;                        % wavelengths [m]

if env.defspectrum == false
    % Parameterized wave height
    if env.defwavesteepness == true
        env.Aw  = env.steepness*(2*pi./env.k)/2; % Wave amplitudes [m]
    end
else
    [env.S_bret,env.moments] = bretschneider(env.Hs,2*pi./env.Tp,env.omega,true);
end

% Calculate group velocity [m/s]
env.Vg = groupvel(env.omega,env.k,env.h);

% Wave TAP per crest width [W/m]
[env] = avgwavepower(env);


% ACE sea states
if solver.calculateACE == 1
    Hs = [2.34; 2.64; 5.36; 2.06; 5.84; 3.26];
    Tp = [7.31; 9.86; 11.52; 12.71; 15.23; 16.50];
    env.ACE.seastates = table(Tp,Hs);
end

for i = 1:nrows
    for j = 1:ncols
        disp([i j])
        % VGOSWEC dimensions:
        w = body.dim.w(i,j);    % current width
        c = body.dim.c(i,j);    % current dist from seafloor
        t = body.dim.t(i,j);    % current thickness
        ht = body.dim.ht(i,j);  % current heights
        
        if body.paramBodyProp == true % If properties are not predefined, use parameterizations:
            % Assign body properties:
            body.prop.I55(i,j) = body.prop.rho_m/3*w*t*ht^3*(1+(t/(2*ht))^2);   % Pitch mass moment of inertia [kg-m^2]
            body.prop.rb(i,j) = ht/2;      % Buoyant torque arm [m]
            body.prop.rg(i,j) = ht/2;      % Gravity torque arm [m]
            body.hydro.C55(i,j) = hs_restoring(body.prop.rho_m,w,t,ht,... % HS Restoring coefficient [Kg-m^2/s^2]
                body.prop.rb(i,j),body.prop.rg(i,j),env.rho,env.g);
            body.prop.mass(i,j) = w*ht*t*body.prop.rho_m;
            body.hydro.Badd(i,j) = 0;
        end
        
        % Hydrodynamic coefficients
        % Michele et al.'s model [1-2]:
        [body.hydro.mu55{i,j},body.hydro.nu55{i,j},body.hydro.X5{i,j},~,...
            body.hydro.mu15{i,j},body.hydro.nu15{i,j},body.hydro.X1{i,j},~,~] ...
            = fun.hydrodynamics.optimizeHydro_M_v2(env.g,env.rho,env.T,  ...
            env.h, w, c, body.prop.I55(i,j), ...
            body.hydro.C55(i,j), solver.n,solver.nmax);
        
    
        % Scale excitation by wave amplitude to obtain forces and moments:
        body.hydro.F5{i,j} = body.hydro.X5{i,j}.*env.Aw;
        body.hydro.F1{i,j} = body.hydro.X1{i,j}.*env.Aw;
        
        % Calculate PTO properties and time-averaged power (TAP) based on
        % control scheme:
        
        % Complex conjugate control:
        if strcmp(body.pto.ctrltype,'CCC')
            if body.pto.pitch_constraint == 1
                % PTO restoring coefficient, CCC case
                body.pto.Cg{i,j} = -(body.hydro.C55(i,j)-env.omega.^2 ...
                    .*(body.prop.I55(i,j)+body.hydro.mu55{i,j}));
                
                % TAP and PTO damping coefficient, CCC case
                [out.TAP{i,j},body.pto.nu_g{i,j}] = CCC_power(env.omega,env.Aw,...
                    body.pto.xi5_max,body.hydro.nu55{i,j},body.hydro.X5{i,j});
            end
            
            % Active damper:
        elseif strcmp(body.pto.ctrltype,'damper')
            % PTO restoring coefficient, damper case
            %             body.pto.Cg{i,j} = body.hydro.C55{i,j}*ones(size(env.omega));
            body.pto.Cg{i,j} = body.pto.Cg_constant*ones(size(env.omega));
            
            % TAPdivA2 and PTO damping coefficient, damper case
            [TAPdivA2,body.pto.nu_g{i,j}] = damper_power(env.rho,env.g,env.Vg,...
                env.k,env.omega,body.prop.I55(i,j),body.hydro.C55(i,j),...
                body.hydro.mu55{i,j},body.hydro.nu55{i,j},body.pto.Cg{i,j});
            
            if env.defspectrum == false
                out.TAP{i,j} = TAPdivA2.*env.Aw.^2; % [W]
            else
                out.TAP{i,j} = abs(trapz(env.omega,TAPdivA2.*2.*env.S_bret)); % [W]
            end
            
            % free
        elseif strcmp(body.pto.ctrltype,'free')
            % Assign zero TAPdivA2 and PTO damping coefficient for free case
            body.pto.nu_g{i,j}  = zeros(size(env.omega));
            body.pto.Cg{i,j}    = zeros(size(env.omega));
            out.TAP{i,j}        = zeros(size(env.omega));
            
        end
        
        % Incorporate simple PTO efficiency
        out.TAP{i,j}= out.TAP{i,j}*body.pto.PTO_eff;
        
        % Capture width ratio
        out.CWR{i,j} = out.TAP{i,j}./(w*env.TAPwave);
        out.maxCWR(i,j) = max(out.CWR{i,j});

        %%% 1-14-21 added Badd to RAO eq
        % Pitch RAO 
        out.RAO{i,j} = (body.hydro.X5{i,j}.*(-(env.omega).^2.*(body.prop.I55(i,j)+body.hydro.mu55{i,j})...
            + 1i*env.omega.*(body.hydro.nu55{i,j}+body.hydro.Badd(i,j)+body.pto.nu_g{i,j}) ...
            + (body.hydro.C55(i,j)+body.pto.Cg{i,j})).^(-1));
        
        % Pitch amplitude
        out.xi5{i,j}      = out.RAO{i,j}.*env.Aw;
        
        
        %%%
        
        env.S_bret
        
        %%%
        
        % Hinge reaction forces (Kurniawan 2012, Tom 2016):
        out.Xr1{i,j} = ((-env.omega.^2.*body.hydro.mu15{i,j}+1i*env.omega.*body.hydro.nu15{i,j})...
            .*out.RAO{i,j} - body.hydro.X1{i,j});
        out.Xr3{i,j} = -(env.rho*w*t*ht - body.prop.mass(i,j))*env.g./env.Aw.*ones(size(env.omega)); % Neglect X3
        out.Fr1{i,j} = out.Xr1{i,j}.*env.Aw;
        out.Fr3{i,j} = out.Xr3{i,j}.*env.Aw;
        out.maxFr1(i,j) = max(out.Fr1{i,j});
        
        if fdn.paramFdnProp == true % If properties are not predefined, use parameterizations:
            % Parameterized foundation dimensions (assumed cylindrical w/ annular cross section):
            fdn.dim.L(i,j) = body.dim.c(i,j); % Length [m]
            out.MbFdnBaseMax(i,j) = max(abs(out.Fr1{i,j}))*fdn.dim.L(i,j);
            
            if fdn.dim.L(i,j) ~= 0
                % fdn.dim.D = sqrt(fdn.prop.Vratio*4/pi*(body.dim.w.*body.dim.ht.*body.dim.t./fdn.dim.L)); % Diameter [m]
                fdn.dim.D(i,j) = (out.MbFdnBaseMax(i,j)/(pi/32*fdn.prop.SigYield/fdn.prop.SF*(1-fdn.prop.gamma^4)))^(1/3);
            else
                fdn.dim.D(i,j) = 0;
            end
            fdn.prop.mass(i,j) = pi/4*(fdn.dim.D(i,j)^2 - (fdn.prop.gamma*fdn.dim.D(i,j))^2)*fdn.dim.L(i,j)*fdn.prop.rho; % Mass [kg]
        end
        
        % Calculate ACE metric:
        if solver.calculateACE == 1
            
            % Calculate ACE time-averaged power by integrating spectrum
            if i==1 && j==1; [out.ACETAP{i,j},env.ACE.S_bret,env.ACE.moments] = calcACETAP(env.omega,TAPdivA2,env.ACE.seastates,1); 
            else;            [out.ACETAP{i,j},~,~] = calcACETAP(env.omega,TAPdivA2,env.ACE.seastates,0); end
            
            % Calculate device and foundation material masses
            steel_mass    = body.prop.mass(i,j)*10^(-3); %[mt]
            concrete_mass = fdn.prop.mass(i,j)*10^(-3);  %[mt]
            out.mtlmass{i,j} = [steel_mass,concrete_mass,0,0,0,0,0];
            
            % Call ACE script (Cole Burge -> NEED TO ADD TO REF!)
            [out.ACCW(i,j),out.CCE{i,j},out.ACE{i,j}] = fun.cost_modeling.calcACE(out.ACETAP{i,j}*10^(-3),out.mtlmass{i,j});
        end
        
    end
end; clear w c t ht concrete_mass steel_mass ACE
toc

%% Time history
%[out,env,body,fdn,param,solver] = fun.hydrodynamics.SSresponse(env,body,fdn,param,solver)


end


%% SUB-FUNCTIONS
% WAVE DISPERSION RELATION SOLVER
function [k] = wave_disp(w,h,g,allowable_error)
k_n = w.^2/g; % make a first guess of k using deep water approximation
err = ones(size(w)); % initialize error

% iterate until max element-wise error is reduced to allowable error:
while max(err) > allowable_error
    w_n = sqrt(g*k_n.*tanh(k_n*h));  % compute new w
    k_n = w.^2./(g*tanh(k_n*h));     % compute new k
    err = abs(w_n - w)./w;           % compute error
end
k = k_n; % assign output as final k
end
%-------------------------------------------------------------------------%
% GROUP VELOCITY
function [Vg] = groupvel(omega,k,h)
Vg = 1/2*omega./k.*(1+2*k*h./(sinh(2*k*h)));
end
%-------------------------------------------------------------------------%
% HYDROSTATIC RESTORING FORCE
function [C55] = hs_restoring(rho_m,w,t,ht,rb,rg,rho,g)
Cbuoy = rho*w*t*ht*rb;
Cgrav = rho_m*w*t*ht*rg;
C55 = (Cbuoy-Cgrav)*g;      
end
%-------------------------------------------------------------------------%
% TIME-AVERAGED WAVE POWER FLUX
function [env] = avgwavepower(env) 
rho = env.rho;
g   = env.g;
k   = env.k;
h   = env.h;
w   = env.omega;

TAPwavecoeff = 1/2*rho*g*(g./k.*tanh(k*h)).^(1/2).*(1+2*k*h./(sinh(2*k*h))); % [W/m^3]

   if env.defspectrum == false
        env.TAPwave = 1/2*env.Aw.^2.*TAPwavecoeff; % [W/m]
   else
        env.TAPwave = abs(trapz(w,TAPwavecoeff.*env.S_bret)); % [W/m]
   end
      
end
%-------------------------------------------------------------------------%
% TAP COEFFICIENT
function [TAP_Coeff] = powercoeff(omega,I55,C55,mu55,nu55,Cg) 

TAP_Coeff = (1 + (1+((C55 + Cg - omega.^2.*(I55 + mu55))./ ...
    (omega.*nu55)).^2).^(1/2)).^(-1); % [-]
end
%-------------------------------------------------------------------------%
% Complex Conjugate Control Scheme Power (REACTIVE CONTROL)
function [TAP,nu_g] = CCC_power(omega,Aw,xi5_max,nu55,X5)

% PTO absorbed power under motion constraints (see [3])
nu_g = ones(length(omega),1);
TAP  = ones(length(omega),1);

Aw = Aw.*ones(length(omega),1);

delta = omega.*xi5_max./Aw * 2.*nu55./X5;

nu_g(delta>=1) = nu55(delta>=1);
nu_g(delta<1) = X5(delta<1)./(omega(delta<1)*xi5_max) - nu55(delta<1);

TAP(delta>=1) = 1/8*Aw(delta>=1).^2.*X5(delta>=1).^2./nu55(delta>=1);
TAP(delta<1) = 1/2*Aw(delta<1).*X5(delta<1).*omega(delta<1)*xi5_max ...
    - 1/2*nu55(delta<1).*omega(delta<1).^2*xi5_max^2;
end
%-------------------------------------------------------------------------%
% UNCONSTRAINED POWER(TIME-VARYING PTO DAMPING, CONSTANT PTO SPRING)
function [TAPdivA2,nu_g] = damper_power(rho,g,Vg,k,omega,I55,C55,mu55,nu55,Cg)
% PTO time averaged power coefficient which ranges b/t 0-0.5
TAP_Coeff = powercoeff(omega,I55,C55,mu55,nu55,Cg);

% PTO damping coefficient
nu_g_coeff = (1+((C55 + Cg - omega.^2.*(I55 + mu55))./(omega.*nu55)).^2).^(1/2);
nu_g = nu55.*nu_g_coeff;

% TAP per wave amplitude squared based on Haskind relation (which relates
% damping and excitation force) [3-4]:
TAPdivA2 = 2*rho*g*Vg./k.*TAP_Coeff; % [W/m^2]
end
%-------------------------------------------------------------------------%
% RAYLEIGH DISTRIBUTION
function [p] = rayleigh(H,variance)
p = H/(4*variance).*exp(-H.^2/(8*variance));
end
%-------------------------------------------------------------------------%
% BRETSCHNEIDER SPECTRUM
function [S,m] = bretschneider(Hs,wm,w,plotout)
Tp = 2*pi/wm;
m = zeros(3,1);
S = 1.25/4*wm^4./(w.^5)*Hs^2.*exp(-1.25*(wm./w).^4);
% S = 5/16*wm.^4./(w.^5)*Hs^2.*exp(-5*wm.^4./(4*w.^4))
m(1) = trapz(w,S);       % zeroth moment m0
m(2) = trapz(w,S.*w);    % first moment m1
m(3) = trapz(w,S.*w.^2); % second moment m2

if plotout == true
    figure
    plot(w,S)
    ylabel('$S_{\eta}(\omega)$ (m\textsuperscript{2}s)','Interpreter','Latex')
    xlabel('$\omega$ (rad/s)','Interpreter','Latex')
end

if S(1) > 0.05*max(S) || ...
        S(end) > 0.05*max(S)
    warning(['The spectral width calculated at 5% of the maximum value of the current seastate ',...
        '(Tp = ', num2str(Tp),' s, Hs = ',num2str(Hs),' m) ',...
        'is outside the range of specified frequencies.',char(10),...
        'Maximum value:',' S(w=w_p=',num2str(2*pi/Tp),' rad/s) = ',num2str(max(S)),' m^2 (5% = ',num2str(0.05*max(S)),' m^2)',char(10),...
        'Value at bound 1: S(w=',num2str(w(1)),' rad/s) = ',num2str(S(1)),' m^2',char(10),...
        'Value at bound 2: S(w=',num2str(w(end)),' rad/s) = ',num2str(S(end)),' m^2',char(10),...
        ])
end


end
%-------------------------------------------------------------------------%
% ACE METRIC TIME-AVG POWER
function [ACETAP,S_bret,moments] = calcACETAP(w,TAPdivA2,seastates,plotspectra)
S_bret = zeros(length(w),height(seastates));
moments = zeros(3,height(seastates));
ACETAP = zeros(1,height(seastates));

for ss = 1:height(seastates)
   
    Hs = seastates.Hs(ss);
    Tp = seastates.Tp(ss);
    wp = 2*pi*Tp^-1;
    % Calculate energy wave spectrum
    [S_bret(:,ss),moments(:,ss)] = bretschneider(Hs,wp,w,false);
    
    % Check bounds at 5% spectral width
    if S_bret(1,ss) > 0.05*max(S_bret(:,ss)) || ...
       S_bret(end,ss) > 0.05*max(S_bret(:,ss))     
        warning(['The spectral width calculated at 5% of the maximum value of seastate '...
            ,num2str(ss),' (Tp = ', num2str(Tp),' s, Hs = ',num2str(Hs),' m) ',...
            'is outside the range of specified frequencies.',char(10),...
            'Maximum value:',' S(w=w_p=',num2str(2*pi/Tp),' rad/s) = ',num2str(max(S_bret(:,ss))),' m^2 (5% = ',num2str(0.05*max(S_bret(:,ss))),' m^2)',char(10),...
            'Value at bound 1: S(w=',num2str(w(1)),' rad/s) = ',num2str(S_bret(1,ss)),' m^2',char(10),...
            'Value at bound 2: S(w=',num2str(w(end)),' rad/s) = ',num2str(S_bret(end,ss)),' m^2',char(10),...
            ])
    end
    
    ACETAP(ss) = abs(trapz(w,TAPdivA2*2.*S_bret(:,ss))); % [W]
end
if plotspectra == 1
    figure()
    plot(w,S_bret)
    legend([repmat('SS ',height(seastates),1),num2str([1:6].')])
end

% Sea state 3 with dir = -70deg set to zero!
ACETAP(3) = 0;
end
%-------------------------------------------------------------------------%
% Longuet-Higgins
% function [joint_p] = LH_jointprob(H,T,moments)
% m0 = moments(1);
% m1 = moments(2);
% m2 = moments(3);
% 
% xi = (H/2)/sqrt(m0)
% Tmean = 2*pi*m0/m1
% nu = sqrt(m2/m0)*Tmean/(2*pi)
% eta = (T-Tmean)/nu
% joint_p = xi.^2/sqrt(2*pi).*exp(-xi.^2.*(1+eta.^2)/2)
% 
% [ximesh,etamesh] = meshgrid(xi,eta);
% 
% surface(ximesh,etamesh,joint_p)
% ylim([0 3])
% end

%% SCRAPS

% [111]:         
%         % Calculate excitation torque
%         body.hydro.X5{i,j} = sqrt(8*env.rho*env.g...
%             *env.Vg.*body.hydro.nu55{i,j}./env.k); % [N]
%         hold on
%         plot(env.omega,body.hydro.X5{i,j})


% omega = 0.1:0.01:1.5;
% T = 2*pi./omega;
% % T = 4:0.1:14;
% Aw = 1;         % wave amplitudes [m]
% Asc = 1;        % wave amplitude scale [-]
% h = 10;         % water depth [m]
% 
% % VGOSWEC dimensions:
% w = 10;         % Width [m]
% t = 0.5;        % Thickness [m]
% c = 0;          % Distance from bottom [m]; c = water depth - device height 
% ht = h-c;
% 
% % Mass properties
% rho_m = rho/2;  % Structural density [kg/m^3]
% rho_m = 36*10^3/72;  % Structural density [kg/m^3]
% I55 = rho_m/3*w*t*ht^3*(1+(t/(2*ht))^2);   % Pitch mass moment of inertia [kg-m^2]
% I55 = 904.4E3;
% rb = ht/2;      % Buoyant torque arm [m]
% rg = ht/2;      % Gravity torque arm [m]
% C55 = (rho*w*t*ht*rb - rho_m*rg)*g;  % Device's buoyant torque calculated from buoyancy [N-m]

% % omega = 0.1:0.01:1.5;
% % T = 2*pi./omega;
% T = 1:0.2:3.8;
% omega = 2*pi./T;
% Aw = [0.004  0.0055	0.0075	0.01	0.0125	0.0155	0.019	0.0225 ...
%       0.0265 0.031	0.0355	0.04	0.0455	0.051	0.0565]'; % wave amplitudes [m]
% Asc = 1;        % wave amplitude scale [-]
% h = 4.5;         % water depth [m]
% 
% % VGOSWEC dimensions:
% w = 0.9398;         % Width [m]
% t = 0.0953;        % Thickness [m]
% c = 4.5 - 0.652;          % Distance from bottom [m]; c = water depth - device height 
% ht = h-c;
% 
% % Mass properties
% rho_m = 25.7/(w*t*ht);  % Structural density [kg/m^3]
% I55 = rho_m/3*w*t*ht^3*(1+(t/(2*ht))^2);   % Pitch mass moment of inertia [kg-m^2]
% I55 = 1.372 + 25.7*(0.651 - 0.364)^2;
% 
% rb = h - 0.294;      % Buoyant torque arm [m]
% rg = h - 0.364;      % Gravity torque arm [m]
% C55 = (rho*w*t*ht*rb - rho_m*rg)*g;  % Device's buoyant torque calculated from buoyancy [N-m]
% % Modified for W2 OSWEC run
% ht = 0.65166875; 
% t = 0.095;
% w = 0.94;
% w_list = w;         % Width [m]
% w2tr = w_list/t;
% c_list = 4.5-ht;          % Distance from bottom [m]; c = water depth - device height

% 
% % Mass properties
% mass = 25.7;    % [kg]
% rho_m = mass/(w_list*t*ht);  % Structural density [kg/m^3]
% I55 = rho_m/3*w*t*ht^3*(1+(t/(2*ht))^2);   % Pitch mass moment of inertia [kg-m^2]
% hinge_d = 0.65166875;
% rg = hinge_d - 0.3624;      % Gravity torque arm [m]
% rb = hinge_d - 0.2940;
% C55 = (rho*w*t*ht*rb - rho_m*rg)*g;  % Device's buoyant torque calculated from buoyancy [N-m]


 
% % Calculate energy wave spectrum
% Hs = 2.34; Tp = 7.31;
% wp = 2*pi*Tp^-1;
% [S_bret,moments] = bretschneider(Hs,wp,env.omega);
% 
% %         4*sqrt(moments(1))
% %         Tz = 2*pi*sqrt(moments(1)/moments(3))
% %         1.408*Tz
% figure
% plot(env.omega,S_bret)
% %
% %         Ai = sqrt(2*S_bret*mean(diff(env.omega)));
% %
% %         trapz(env.omega,S_bret)
% %         plot(env.omega,Ai)
% 
% body.pto.Cg{i,j} = body.hydro.C55{i,j}*ones(size(env.omega));
% [TAPdivA2,body.pto.nu_g{i,j}] = unconstr_power_spec(env.rho,env.g,env.Vg,...
%     env.k,env.omega,body.prop.I55(i,j),body.hydro.C55{i,j},...
%     body.hydro.mu55{i,j},body.hydro.nu55{i,j},body.pto.Cg{i,j});
% 
% TAP1 = TAPdivA2 * (Hs/2).^2; % [W]
% TAP2 = trapz(env.omega,TAPdivA2*2.*S_bret); % [W]
% plot(env.omega,TAP1); hold on
% yline(TAP2)
% 
% %             % Calculate wave height distribution
% %
% %             WaveHeights = 0:0.05:2.3*Hs
% %             variance = (Hs/4)^2
% %             p = rayleigh(WaveHeights,variance)
% %             plot(WaveHeights,p)
% %
% %             nn = TAPdivA2*WaveHeights
% %
% %             TAPw = []
% %             for row = 1:size(nn,1)
% %                 %            TAPw(row,1) = trapz(env.omega(row),nn(row,:).*p)
% %                 TAPw2(row,1) = sum(nn(row,:).*p*mean(diff(env.omega)));
% %
% %
% %             end
% %             figure
% %             plot(env.omega,TAPw2); hold on
% %             yline(TAP2)
% 
% 
% 
% %         f = env.omega/(2*pi)
% %         T = f.^-1
% %
% %         joint_p = LH_jointprob(WaveHeights,(env.omega/(2*pi)).^(-1),moments)

%-------------------------------------------------------------------------%
% function optimizeHydro_RD(parallel, g, rho, Tp, A1p, Ap, hp, wp, cp, Ip, Cp, n, np)
%   
%	Input: 1) Environmental:        - T (wave periods)
%                                   - A1 (wave amplitude)
%                                   - h Water depth
%                                   - A Wave scale amplitude

%             Geometric properties: - w Device's width
%                                   - c Distance from seafloor
%                                   - I Device's moment of inertia
%                                   - C Device's Buoyancy Torque
%
%             Others:               - parallel on/off
%                                   - g gravitational constant
%                                   - rho water density
%                                   - n and nmax: numerical convergence
%                                   criteria.
%
%   Output: - miu (added inertia)
%           - nu (radiation damping)
%           - F (exciting force)
%           - Capture Factor
%           - Power
%           - nu_ptop optimized PTO damping
%-------------------------------------------------------------------------%
%   Based on Renzi and Dias (2013) - Hydrodynamics of the oscillating wave
%                                    surge converter in the open ocean.
%
%   Written by: Jessica Nguyen, PhD 
%               University of Massachusetts Amherst
%               nvnguyen@umass.edu
%-------------------------------------------------------------------------%
function [miu, nu, nu_pto, F, P, Cf_opt, Cg] = optimizeHydro_RD_v2(parallel,...
                                                    g, rho, T, Aw, Asc, h,...
                                                    w, c, I, C, n, nmax)
                                                
%% Setup for constants and other parameters (NO Modification Needed); 
% Environmental variables;
omega = 2*pi./T;
Aw = Aw/Asc;

% Setup empty holders;
ncases = length(omega);
k0 = zeros(ncases, 1);
kn = zeros(ncases,n);
knTemp = zeros(ncases, 2*(n-1));

f0 = zeros(ncases,1);
d0 = zeros(ncases,1);
fn = zeros(ncases,n);
fnTemp = zeros(ncases,2*(n-1));

Up = zeros(nmax,nmax);
alpha = zeros(n, nmax);
alpha0n = zeros(n, ncases);
miu = zeros(1,ncases); 
nu = zeros(1,ncases); 
nu_pto = zeros(1,ncases); 
F = zeros(1,ncases); 

%% Calculates wavenumbers;
for i = 1:ncases 
	err = 1; ktemp = 1;
    while (err >= 1e-12)
        func = omega(i)^2 - g*ktemp*tanh(ktemp*h);
        dfunc = -g*(tanh(ktemp*h) + h*ktemp*(sech(ktemp*h))^2);
        ktemp_new = ktemp - func/dfunc; 
        err = abs(ktemp_new-ktemp);
        ktemp = ktemp_new;
    end  
    k0(i) = ktemp;
    [f0(i), d0(i)] = cal_fdn(k0(i), h, c, g, omega(i));
    
    for j = 1:2:2*(n-1)
        func = @(x) (omega(i)^2 + g*x*tan(x*h));
        knTemp(i,j) = -1i*fzero(func, [(j*pi/2+j*7.5e-4)/h (j+1)*pi/2/h]);
        [fnTemp(i,j), ~] = cal_fdn(knTemp(i,j), h, c, g, omega(i));
    end
end

for i = 1:ncases
    kn(i,:) = [k0(i) nonzeros(knTemp(i,:))'];
    fn(i,:) = [f0(i) nonzeros(fnTemp(i,:))'];
end

%% Calculates Chebyshev coefficients;
% Calculates vo (zeros of the first kind Chebyshev polynomials);
j = 0:1:nmax-1;
vo = cos((2*j+1)*pi/2/(nmax+1));

% Calculates Chebyshev functions;
Up(:,1) = 1;
Up(:,2) = 2*vo;
for i = 1:nmax
    for j = 3:nmax
        dn = Up(i,j-1); dnm1 = Up(i,j-2);
        Up(i,j) = 2*vo(i)*dn - dnm1; 
    end
end 

% Form Cpn matrix and solve for alpha coefficients for each wave period;   
if (parallel == 1)      
    parfor i = 1:ncases       % Lopps through each wave period;
        disp(i);
        alpha = zeros(n, nmax);
        knRow = kn(i,:);
        fnRow = fn(i,:);
        for j = 1:n      % Loops through each kn (or fn);
            alpha(j, :) = cal_Alpha_Beta(nmax, vo, Up, knRow(j), fnRow(j), w);
        end
        alpha0n(:,i) = alpha(:,1);
    end
else
    for i = 1:ncases          % Lopps through each wave period;
        disp(i);
        for j = 1:n      % Loops through each kn (or fn);
            alpha(j, :) = cal_Alpha_Beta(nmax, vo, Up, kn(i,j), fn(i,j), w);
        end
        alpha0n(:,i) = alpha(:,1);
    end
end

%% Computes added inertia, radiation damping, and exciting torque;
% Evaluates the normalized added inertia and damping coefficients;
for i = 1:ncases
    miu(i) = rho*w*pi/4*real(sum(fn(i,:).*alpha0n(:,i)')); 
    nu(i) = rho*w*omega(i)*pi/4*fn(i,1)*imag(alpha0n(1,i));
    F(i) = rho*w*1i*omega(i)*pi/4*Aw*alpha0n(1,i)*d0(i);
end    

%% Computes average power and capture factor;
for i = 1:ncases
    nu_pto(i) = sqrt((C - (I + miu(i))*omega(i)*omega(i))^2/...
                    (omega(i)*omega(i)) + nu(i));    
end

Cg = omega'./2./k0.*(1 + 2*k0*h./sinh(2*k0*h)); 
P = abs(F).*abs(F)./4./(nu + nu_pto);
Cf_opt = P./(1/2)/Aw/Aw./Cg; 

% Maximum power
% nu_ptop_opt = nup;
% Pmax = abs(Fp).*abs(Fp)./4./(nup + nu_ptop_opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------SUB-FUNCTIONS----------------------------%

%% Calculates fn and dn;
function [f, d] = cal_fdn(k, h, c, g, omega)

f = (sqrt(2)*(k*(h-c)*sinh(k*h) + cosh(k*c) - cosh(k*h)))/...
        (k*k*sqrt(h + g*(sinh(k*h))^2/omega/omega));
d = g*k*sqrt(h + g*(sinh(k*h))^2/omega/omega)/sqrt(2)/omega/cosh(k*h);

%% Calculates the coefficients for Chebyshev second kinds
function [alpha, Cpn] = cal_Alpha_Beta(np, voVec, Up, kn, fn, w)

% For wave fn (n = 0,1,2,...,n), forms the Cpn np x np matrix;
Cpn = zeros(np,np);
for i = 1:np % Horizontal loop calculate the coefficients of the Cpn matrix;
    vo = voVec(i);
    Rn = form_Rn(np, kn, vo, w);
    
    % Vertical loop calculate the coefficients of the Cpn matrix;
    for j = 1:np
        term1 = -pi*(j)*Up(i,j);       % Alternatively: -(j)*chebyshevU(j-1,vo);
       
        if (j == 2)
            term1 = 0.5*term1;
        end
        
        fun = @(u) sqrt(1-u.^2).*chebyshevU(j-1,u).*Rn(u)./abs(vo - u);%*chebyshevU(j-1,u).
        term2 = 1i*kn*pi*w/4*integral(fun, -1, 1);

        Cpn(i,j) = term1 + term2; 
    end
end 

% Solves for the coefficients alpha_pn values; 
RHS = pi*w*fn*ones(np, 1);
alpha = -Cpn\RHS;

%% Calculate the remainders in Hankel function
function Rn = form_Rn(np, kn, vo, w)

% Constructs the first and second terms;
t1 = @(u) besselj(1, kn*w/2*abs(vo - u));
t2 = @(u) 1+ 2*1i/pi*(log(kn*w/4*abs(vo-u)) + 0.5772156649015329);

% Calculates the last term;
sumJ = @(u) 0;
for j = 2:2*np
    q = 1:j-1;
    sumQ = sum(2./q);
    cs = (-1)^(j+1)/factorial(j)/factorial(j-1)*(1/j + sumQ);
    ks = @(u) (kn*w/4*abs(vo - u)).^(2*j-1);
    sumJ = @(u) sumJ(u) + cs*ks(u);
end

t3 = @(u) 1i/pi*(kn*w/4*abs(vo-u) + sumJ(u));

% Forms the remainder Rn equation;
Rn = @(u) t1(u).*t2(u) - t3(u);









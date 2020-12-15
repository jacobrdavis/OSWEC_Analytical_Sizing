%-------------------------------------------------------------------------%
% function optimizeHydro_M(g, rho, T, Aw, Asc, h, w, c, I, C, n, nmax)
%   
%	Input: 1) Environmental:        - T (wave periods)
%                                   - Aw (wave amplitude)
%                                   - h Water depth
%                                   - Asc Wave scale amplitude

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
%   Based on Michele et al. (2016) - Theory of the synchronous motion of an 
%                                    array of ﬂoating ﬂap gates oscillating
%                                    wave surge converter.
%
%   Written by: Jessica Nguyen, PhD 
%               University of Massachusetts Amherst
%               nvnguyen@umass.edu
%-------------------------------------------------------------------------%
function [mu55, nu55, X5, P55, mu15, nu15, X1, nu_pto55, Cf_opt55, Cg55] = ...
                optimizeHydro_M_v2(g, rho, T, Aw, Asc, h, w, c, I, C, n, nmax) 

%% Setup for constants and other parameters (NO Modification Needed); 
% Environmental variables;
omega = 2*pi./T;
Aw = Aw/Asc;

% Radial Mathieu orders and parameters;
Xi = 0;      
order = 1:2:nmax;

% Setup empty holders;
ncases = length(omega);
k0 = zeros(ncases, 1);
kn = zeros(ncases,n);
knTemp = zeros(ncases, 2*(n-1));

d0 = zeros(ncases,1);
f0 = zeros(ncases,1);
fn = zeros(ncases,n);
fnTemp = zeros(ncases,2*(n-1));

lambda0 = zeros(ncases,1);
lambdan = zeros(ncases,n);
lambdaTemp = zeros(ncases,2*(n-1));

mu55 = zeros(ncases, 1);
nu55 = zeros(ncases, 1);
X5 = zeros(ncases, 1);
nu_pto55 = zeros(ncases, 1);

mu15 = zeros(ncases, 1);
nu15 = zeros(ncases, 1);
X1 = zeros(ncases, 1);

%% Calculates wavenumbers and solutions of the dispersion relation;
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
    lambda0(i) = cal_lambdan(k0(i), h, c, g, omega(i));
    
    for j = 1:2:2*(n-1)
        func = @(x) (omega(i)^2 + g*x*tan(x*h));
        knTemp(i,j) = -1i*fzero(func, [(j*pi/2+j*7.5e-4)/h (j+1)*pi/2/h]);
        [fnTemp(i,j), ~] = cal_fdn(knTemp(i,j), h, c, g, omega(i));
        lambdaTemp(i,j) = cal_lambdan(knTemp(i,j), h, c, g, omega(i));
    end
end

for i = 1:ncases
    kn(i,:) = [k0(i) nonzeros(knTemp(i,:))'];
    fn(i,:) = [f0(i) nonzeros(fnTemp(i,:))'];
    lambdan(i,:) = [lambda0(i) nonzeros(lambdaTemp(i,:))'];
end
tau_n = (w*w*kn.*kn/16);

%% Calculates added inertia, damping coefficients, and exciting force;
for i = 1:ncases
    miu_temp55 = 0; miu_temp15 = 0;
    for j = 1:n
        q = tau_n(i,j);
        [~, AmBm,~] = eig_AmBm(nmax, q);

        vtemp = 0;
        for k = 1:length(order)
            Apm = AmBm(:,(order(k)+1)/2);  
            npm =  Npm(nmax, Xi, q, Apm);
            dhpm = dHpm(nmax, Xi, q, Apm);
            
            vtemp = vtemp + Apm(1)*Apm(1)*npm/4./dhpm;
        end

        miu_temp55 = miu_temp55 + fn(i,j)*fn(i,j)*imag(vtemp);
        miu_temp15 = miu_temp15 + lambdan(i,j)*fn(i,j)*imag(vtemp);
        if (j == 1)
            nu_temp = real(vtemp);
            f_temp = vtemp;
        end
    end
    mu55(i) = rho*w*w*pi*miu_temp55;
    nu55(i) = -rho*omega(i)*w*w*f0(i)*f0(i)*pi*nu_temp;
    X5(i) = abs(rho*omega(i)*Aw*w*w*f0(i)*d0(i)*pi*f_temp);
    
    mu15(i) = rho*w*w*pi*miu_temp15;
    nu15(i) = -rho*omega(i)*w*w*lambda0(i)*f0(i)*pi*nu_temp;
    X1(i) = abs(rho*omega(i)*Aw*w*w*lambda0(i)*d0(i)*pi*f_temp);    
end

%% Computes average power and capture factor;
for i = 1:ncases
    nu_pto55(i) = sqrt((C - (I + mu55(i))*omega(i)*omega(i))^2/...
                    (omega(i)*omega(i)) + nu55(i));    
end

Cg55 = omega'./2./k0.*(1 + 2*k0*h./sinh(2*k0*h)); 
P55 = abs(X5).*abs(X5)./4./(nu55 + nu_pto55);
Cf_opt55 = P55./(1/2)/Aw/Aw./Cg55; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------SUB-FUNCTIONS----------------------------%  
%% Calculates fn and dn;
function [f, d] = cal_fdn(k, h, c, g, omega)
f = (sqrt(2)*(k*(h-c)*sinh(k*h) + cosh(k*c) - cosh(k*h)))/...
        (k*k*sqrt(h + g*(sinh(k*h))^2/omega/omega));
d = g*k*sqrt(h + g*(sinh(k*h))^2/omega/omega)/sqrt(2)/omega/cosh(k*h);

%% Calculates lambda;
function lambda = cal_lambdan(k, h, c, g, omega)
lambda = (sqrt(2)*(sinh(k*h)-sinh(k*c)))/...
        (k*sqrt(h + g*(sinh(k*h))^2/omega/omega));

%% Calculates coefficients in Mathieu Functions Type 4: odd-odd;
function [alpha, AmBm, order] = eig_AmBm(nmax, q)
% Constructs the LHS matrix;
nk = 0:nmax-1;
nonDiag = q*ones(1,length(nk)-1);
M = diag((2*nk+1).^2) + diag(nonDiag, -1) + diag(nonDiag, 1);
M(1,1) = M(1,1)-q;

% Finds the eigenvalues and eigenvectors;
[eigVec, eigvalues] = eig(M,'nobalance');
[alpha, num] = sort(diag(eigvalues));     % alpha is the eigenvalue vector;
eigV = eigVec(:,num);                     % matrix of eigenvectors
order = 2*nk+1;                           % odd values of order n

%   Compute matrix mv of processed eigenvectors
AmBm = zeros(nmax, size(eigV,2));
for k = 1:size(eigV,2)
    eigVcol = eigV(:,k).*order';         
    norm_eigVcol = eigV(:,k)/sum(eigVcol);       
    scale = sum(norm_eigVcol.*norm_eigVcol);
    AmBm(:,k) = norm_eigVcol./sqrt(scale);          
end

%% Evaluates angular Mathieu solutions and their derivatives for Type 4: odd-odd;;
function ceSe = ceSepm(nmax, eta, Apm)
ceSe = 0;
for j = 1:nmax
    ceSe = ceSe + Apm(j)*sin((2*j-1)*eta);
end

function dceSe = dceSepm(nmax, eta, Apm)
dceSe = 0;
for j = 1:nmax
    dceSe = dceSe + (2*j-1)*Apm(j)*cos((2*j-1)*eta);
end

%% Evaluates Radial Mathieu solutions and their derivatives for Type 4: odd-odd;
function npm = Npm(nmax, Xi, q, Apm)
v1 = sqrt(q)*exp(-Xi); 
v2 = sqrt(q)*exp(Xi);

npm = 0;
for j = 0:nmax-1
    npm = npm + (-1)^j*Apm(j+1)*(besselj(j,v1)*bessely(j+1,v2) - besselj(j+1,v1)*bessely(j,v2));
end

so = ceSepm(nmax, pi/2, Apm)*dceSepm(nmax, 0, Apm);
npm = so/sqrt(q)/Apm(1)/Apm(1)*npm;

function dnpm = dNpm(nmax, Xi, q, Apm)
v1 = sqrt(q)*exp(-Xi); 
v2 = sqrt(q)*exp(Xi);

dnpm = 0;    
for j = 0:nmax-1
    dnpm = dnpm + (-1)^j*Apm(j+1)*((v2 + v1)*(besselj(j,v1)*bessely(j,v2) + besselj(j+1,v1)*bessely(j+1,v2))...
        -(2*j+1)*(besselj(j+1,v1)*bessely(j,v2) + besselj(j,v1)*bessely(j+1,v2)));
end

so = ceSepm(nmax, pi/2, Apm)*dceSepm(nmax, 0, Apm);
dnpm = so/sqrt(q)/Apm(1)/Apm(1)*dnpm;


function djpm = dJpm(nmax, Xi, q, Apm)
v1 = sqrt(q)*exp(-Xi); 
v2 = sqrt(q)*exp(Xi);
u = v2-v1;

djpm = 0;
for j = 0:nmax-1
    djpm = djpm + Apm(j+1)*(besselj(2*j,u) - besselj(2*j+2,u));
end

djpm = dceSepm(nmax, 0, Apm)/Apm(1)*cosh(Xi)*djpm;

function dhpm = dHpm(nmax, Xi, q, Apm)

djpm = dJpm(nmax, Xi, q, Apm);
dnpm = dNpm(nmax, Xi, q, Apm);

dhpm = djpm + 1i*dnpm;



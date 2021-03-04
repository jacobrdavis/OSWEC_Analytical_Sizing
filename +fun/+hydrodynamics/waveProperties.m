function [projAreaKEperSA] = waveProperties(env,body,i,j)


% phi = 0 or pi/2;
% A_complex = Aw*exp(1i*phi)
% function [eta] = eta_function(A_complex,k,omega,x,t)
% eta = A_complex*exp(1i*(k*x - omega*t));
% end

env.rho;
h = env.h;


% Do this later w/o loop!
for w = 1:length(env.omega)
L = 2*pi./env.k(w);
k = env.k(w);

Aw    = env.Aw(w);
omega = env.omega(w);
T = 2*pi/omega;

% nx = 2000;
% dx = L/nx;
x = 0;
% x2 = L;
% x = 0:dx:L;

t2 = T;
t1 = 0;

% nz = 2000;
z2 = body.dim.ht(i,j) + body.dim.c(i,j) - env.h;
z1 = body.dim.c(i,j) - env.h; 
% dz = (z2-z1)/nz;
% z = transpose(z1:dz:z2);


% [x_grid, z_grid] = meshgrid(x,z);


% [Dean and Dalrymple p.97, eq. 4.73]

% numerical validation
% F = 1/2*(cosh(2*k*(env.h+z_grid))+cos(2*(k*x_grid-omega*t)));
% avgKE(w) = C*trapz(x,trapz(z,F,2));
C = env.rho/(2)*(env.g*Aw*k/omega*1/cosh(k*env.h))^2;
F = 1/2*(t2-t1)^(-1)*((t2-t1)/(2*k)*(sinh(2*k*(h+z2)) - sinh(2*k*(h+z1)))...
            -(z2-z1)/(2*omega)*(sin(2*(k*x-omega*t1)) - sin(2*(k*x-omega*t2))));

projAreaKEperSA(w,1) = C*F;   %[J/m^2]
projAreaKEperW(w,1) = C*F*L; %[J/m]
fullKEperSA(w,1) = 1/16*env.rho*env.g*(2*Aw)^2;
fullKEperW(w,1) = fullKEperSA(w)*L;


% if omega== env.omega(abs(env.omega-7) == min(abs(env.omega-7)))
% Fcheck = -1/2*env.rho*env.g*(z2^2-z1^2) + env.rho*env.g*Aw/k*cos(k*x1-omega*t)/cosh(k*h)*(sinh(k*(h+z2))-sinh(k*(h+z1)))
% z = linspace(z1,z2,100)
% 
% Fchecknum = trapz(z,(-env.rho*env.g*z + env.rho*env.g*Aw*cosh(k*(h+z))/cosh(k*h)))*body.dim.w(i,j).*cos(k*x1-w*t)
%end
end



%     function [xi,xi_dot,zeta,zeta_dot] = particle_trajectories(A,k,omega,h,z,x,t)
%         
%         
%         
%         xi_amp = -A.*cosh(k.*(z+h))./sinh(k.*h);
%         xi = xi_amp.*sin(k.*x-omega.*t);
%         xi_dot = -omega*xi_amp.*cos(k.*x-omega.*t + pi/2);
%         
%         
% %         xi_amp(1,2,1)*sin(k*x(1,2,1)-omega*t(1,2,1))
%         
%         zeta_amp = A.*sinh(k.*(z+h))./sinh(k.*h);
%         zeta = h+z + zeta_amp.*cos(k.*x-omega.*t + pi/2);
%         zeta_dot = omega*zeta_amp.*sin(k.*x-omega.*t + pi/2);
%             
%     end



end


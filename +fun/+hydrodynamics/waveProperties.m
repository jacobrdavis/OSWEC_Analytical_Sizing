function [projAreaKEperW] = waveProperties(env,body,i,j)


% phi = 0 or pi/2;
% A_complex = Aw*exp(1i*phi)
% function [eta] = eta_function(A_complex,k,omega,x,t)
% eta = A_complex*exp(1i*(k*x - omega*t));
% end

env.rho;
h = env.h;
t = 0;

% Do this later w/o loop!
for w = 1:length(env.omega)
L = 2*pi./env.k(w);
k = env.k(w);

Aw    = env.Aw(w);
omega = env.omega(w);


% nx = 2000;
% dx = L/nx;
x1 = 0;
x2 = L;
% x = 0:dx:L;

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
C = env.rho/(2*L)*(env.g*Aw*k/omega*1/cosh(k*env.h))^2;
F = 1/2*((x2-x1)/(2*k)*(sinh(2*k*(h+z2)) - sinh(2*k*(h+z1)))...
            +(z2-z1)/(2*k)*(sin(2*(k*x2-omega*t)) - sin(2*(k*x1-omega*t))));

projAreaKEperSA(w,1) = C*F;   %[J/m^2]
projAreaKEperW(w,1) = C*F*L; %[J/m]
fullKEperSA(w,1) = 1/16*env.rho*env.g*(2*Aw)^2;
fullKEperW(w,1) = fullKEperSA(w)*L;
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


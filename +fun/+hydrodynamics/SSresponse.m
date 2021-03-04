function [out,env,body,fdn,param,solver] = SSresponse(env,body,fdn,param,solver)
i=1;j=1;
w = env.omega(10)
F0 = env.Aw(env.omega == w)*body.hydro.X5{i,j}(env.omega == w);
M = body.prop.I55(i,j) + body.hydro.mu55{i,j}(env.omega == w);
B = body.hydro.Badd(i,j) + body.hydro.nu55{i,j}(env.omega == w) + body.pto.nu_g{i,j}(env.omega == w);
C = body.hydro.C55(i,j) + body.pto.Cg{i,j}(env.omega == w);
fs = 100;
y0 = [0; 0];
tspan = 0:fs^-1:100-fs^-1;
[t,y] = ode45(@(t,y)my_ode(t,y,F0,w,M,B,C),tspan,y0);
 unitstep = t<=t(end)/2;
ramp = t.*unitstep
plot(t,ramp)
x = y(:,1);
plot(t,abs(x)); hold on
[pks,loc] = findpeaks(y(:,1));
plot(t(loc(2:end)),pks(2:end),'o')

xlabel('Time'); ylabel('x')
xlim([min(t) max(t)])

T_est = mean(diff(t(loc)));
omega_est = 2*pi*T_est^-1; 


function dydt = my_ode(t,y,F0,w,M,B,C)
dydt = [
   y(2);
   (real(F0*exp(1i*(w*t)))-B*y(2)-C*y(1)).*M^(-1)];
end



end


%---------------------------- Header -------------------------------------%
% Title: VGOSWEC Sizing Tool W2 Post Processing Script
% Author: J. Davis
% Date: 12-15-20
% Last Update: 12-15-20
%-------------------------------------------------------------------------%
%% Settings
close all

% Indices of result to plot:
i=1;j=1;

% Define normalizations and labels:
omega_normalization = 1;    omega_label = '$\omega$ (rad s\textsuperscript{-1})';
% omega_normalization = sqrt(env.h/env.g); omega_label = '$\omega* = \omega (h/g)^{1/2} [-]$'; % [rad/s]/[rad/s]
mu_normalization    = 1;	mu_label    = '$\mu_{55} [Kgm^2]$';
nu_normalization    = 1;	nu_label    = '$\nu_{55} [Kgm^2/s]$';
X5_ma_normalization = 1;	X5_ma_label = 'Ex Torque Magnitude $|X_5| [Nm]$';
X5_ph_normalization = 1;	X5_ph_label = 'Ex Torque Phase $\varphi_5 [rad]$';
X1_ma_normalization = 1;	X1_ma_label = 'Ex Force Magnitude $|X_1| [N]$';
X1_ph_normalization = 1;	X1_ph_label = 'Ex Force Phase $\varphi_1 [rad]$';
C55_normalization   = 1;	C55_label   = 'Linear Restoring Coefficient $C_{55} [Kgm^2/s^2]$';
xi5_normalization   = 1;	xi5_label   = 'Pitch Amplitude $\xi_{5} [deg]$';
RAO_normalization   = 1;	RAO_label   = 'RAO $\xi_{5}/A$';
F1_normalization    = 1;	F1_label    = '$|F_1| [N]$';
Fr1_normalization   = 1;	Fr1_label   = 'Hinge Surge Reaction Force Magnitude $[N]$';


%% Plots

% Added mass and radiation damping
figure
colororder({'k','#50F'})
yyaxis left;plot(env.omega*omega_normalization, ...
    body.hydro.mu55{i,j}*mu_normalization,'LineWidth',1.25,...
    'DisplayName','Analytical Added Mass');hold on
ylabel(mu_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
plot(Val.Data1.w,squeeze(Val.Data1.A(5,5,:))*Val.Data1.rho,...
    '--','Color','k','DisplayName',[Val.Source1,' Added Mass'])
yyaxis right;
plot(2*pi./env.T, body.hydro.nu55{i,j}.*nu_normalization,'LineWidth',1.25,...
    'DisplayName','Analytical Radiation Damping');
ylabel(nu_label,'Interpreter','Latex');
plot(Val.Data1.w,squeeze(Val.Data1.B(5,5,:))*Val.Data1.rho.*Val.Data1.w.',...
    '--','Color','#50F','DisplayName',[Val.Source1,' Radiation Damping'])
legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',16)

% Linear restoring coefficient
figure
colororder({'k','#50F'})
plot(env.omega*omega_normalization, ...
    body.hydro.C55*ones(size(env.omega)),'LineWidth',1.25,...
    'DisplayName','Analytical Linear Restoring Coefficient');hold on
ylabel(C55_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
ylim([0 round(1.1*max([body.hydro.C55  Val.Data1.C(5,5)*Val.Data1.rho*Val.Data1.g]))])
plot(Val.Data1.w,Val.Data1.C(5,5)*ones(size(Val.Data1.w))*Val.Data1.rho*Val.Data1.g,...
    '--','Color','k','DisplayName',[Val.Source1,' Linear Restoring Coefficient'])
legend('Interpreter','Latex','Location','SouthEast')
set(gca,'FontSize',16)

% Excitation torque
% Magnitude
figure
plot(env.omega*omega_normalization,abs(body.hydro.X5{1,1})*X5_ma_normalization,'k',...
    'DisplayName','Analytical Excitation Moment Magnitude'); hold on
plot(Val.Data1.w,squeeze(Val.Data1.ex_ma(5,1,:))*Val.Data1.rho.*Val.Data1.g,...
    '--','Color','k','DisplayName',[Val.Source1,' Excitation Moment Magnitude'])
xlabel(omega_label,'Interpreter','Latex')
ylabel(X5_ma_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
legend('Interpreter','Latex','Location','SouthEast')
set(gca,'FontSize',16)
% Phase
figure
plot(env.omega*omega_normalization,pi-angle(body.hydro.X5{1,1})*X5_ph_normalization,'k',...
    'DisplayName','Analytical'); hold on
plot(Val.Data1.w,squeeze(Val.Data1.ex_ph(5,1,:)),...
    '--','Color','k','DisplayName',[Val.Source1])
xlabel(omega_label,'Interpreter','Latex')
ylabel(X5_ph_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
legend('Interpreter','Latex','Location','SouthEast')
set(gca,'FontSize',16)

% Excitation force
% Magnitude
figure
plot(env.omega*omega_normalization,abs(body.hydro.X1{1,1})*X1_ma_normalization,'k',...
    'DisplayName','Analytical Excitation Force'); hold on
plot(Val.Data1.w,squeeze(Val.Data1.ex_ma(1,1,:))*Val.Data1.rho.*Val.Data1.g,...
    '--','Color','k','DisplayName',[Val.Source1,' Excitation Force'])
xlabel(omega_label,'Interpreter','Latex')
ylabel(X1_ma_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
legend('Interpreter','Latex','Location','SouthEast')
set(gca,'FontSize',16)
% Phase
figure
plot(env.omega*omega_normalization,pi-angle(body.hydro.X1{1,1})*X1_ph_normalization,'k',...
    'DisplayName','Analytical'); hold on
plot(Val.Data1.w,squeeze(Val.Data1.ex_ph(1,1,:)),...
    '--','Color','k','DisplayName',[Val.Source1])
xlabel(omega_label,'Interpreter','Latex')
ylabel(X1_ph_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
legend('Interpreter','Latex','Location','SouthEast')
set(gca,'FontSize',16)

% Pitch amp xi5
figure
plot(env.omega,abs(out.xi5{i,j})*180/pi,'DisplayName','Analytical'); hold on
plot(Val.Data2.w,Val.Data2.xi5Amp*180/pi,'-',...
    'DisplayName',[strrep(Val.Source2,'_',' ')])    
plot(Val.Data3.w,Val.Data3.phi,'^',...
    'DisplayName',[strrep(Val.Source3,'_',' ')]) 
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
ylabel(xi5_label,'Interpreter','Latex')
set(gca,'FontSize',16)
legend('Interpreter','Latex')

% 
% 
% plot(Val.Data1.w,squeeze(Val.Data1.ex_ph(5,1,:)))
% exph = squeeze(Val.Data1.ex_ph(5,1,:))
% 
% figure
% yyaxis left
% plot(Val.Data1.w,squeeze(Val.Data1.ex_ma(5,1,:))); hold on
% plot(Val.Data1.w,squeeze(Val.Data1.ex_ma(5,1,:)).*exp(1i*squeeze(Val.Data1.ex_ph(5,1,:)))); hold on
% yyaxis right
% plot(Val.Data1.w,squeeze(Val.Data1.ex_ph(5,1,:)))


% RAO
figure
plot(env.omega,abs(out.RAO{1,1}),'DisplayName','Analytical'); hold on
% plot(env.omega,RAO_wph,'DisplayName','Analytical'); hold on
plot(Val.Data2.w,Val.Data2.xi5Amp./(Val.Data2.H/2),'-',...
    'DisplayName',[Val.Source2,' RAO'])
plot(Val.Data3.w,Val.Data3.RAO.*Val.Data3.k,'s',...
    'DisplayName',[Val.Source3,' RAO'])
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
ylabel(RAO_label,'Interpreter','Latex')
set(gca,'FontSize',16)
legend('Interpreter','Latex')

% Surge force
figure
plot(env.omega*omega_normalization,abs(body.hydro.F1{i,j}),'k',...
    'DisplayName','Analytical'); hold on
plot(Val.Data2.w,Val.Data2.Fr1Amp,'-',...
    'DisplayName',[Val.Source2])
xlabel(omega_label,'Interpreter','Latex')
ylabel(F1_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
set(gca,'FontSize',16)
legend('Interpreter','Latex')

% Hinge reaction forces
figure
plot(env.omega,abs(out.Fr1{i,j}),'DisplayName','Analytical Hinge Reaction Force'); hold on
plot(Val.Data2.w,Val.Data2.Fcnstr1Amp,'-',...
    'DisplayName',[Val.Source2,' Hinge Reaction Force'])
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
ylabel(Fr1_label,'Interpreter','Latex')
set(gca,'FontSize',16)
legend('Interpreter','Latex')



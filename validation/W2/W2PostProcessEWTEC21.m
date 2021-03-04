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
mu55_normalization    = 1;	mu55_label    = '$\mu_{55}$ (kg-m\textsuperscript{2})';
mu15_normalization    = 1;	mu15_label    = '$\mu_{15}$ (kg-m)';
nu55_normalization    = 1;	nu55_label    = '$\nu_{55}$ (kg-m\textsuperscript{2}s\textsuperscript{-1})';
nu15_normalization    = 1;	nu15_label    = '$\nu_{15}$ (kg-m s\textsuperscript{-1})';
X5_ma_normalization = 1;	X5_ma_label = '$|X_5|$ (kN-m)';
X5_ph_normalization = 1;	X5_ph_label = '$\varphi_5$ (rad)';
X1_ma_normalization = 1;	X1_ma_label = '$|X_1|$ (kN)';
X1_ph_normalization = 1;	X1_ph_label = '$\varphi_1$ (rad)';
C55_normalization   = 1;	C55_label   = 'Linear Restoring Coefficient $C_{55} [Kgm^2/s^2]$';
xi5_normalization   = 1;	xi5_label   = 'Pitch Amplitude $\xi_{5} [deg]$';
RAO_normalization   = 1;	RAO_label   = 'RAO $\xi_{5}/A$';
F1_normalization    = 1;	F1_label    = '$|F_1| [N]$';
Fr1_normalization   = 1;	Fr1_label   = 'Hinge Surge Reaction Force Magnitude $[N]$';

omega_limits = [min(env.omega) max(env.omega)]*omega_normalization;
% omega_limits = [min(env.omega) 10]*omega_normalization;

%% Plots

% Added mass 55
figure()
x0=4; y0=4;
width=3.5; height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
plot(env.omega*omega_normalization, ...
    body.hydro.mu55{i,j}*mu55_normalization,'-','LineWidth',1.5,...
    'Color','k','DisplayName','Analytical');hold on
plot(Val.Data1.w,squeeze(Val.Data1.A(5,5,:))*Val.Data1.rho,...
    '--','Color','k','DisplayName','WAMIT')
ylabel(mu55_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim(omega_limits)
%legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',12)

% radiation damping 55
figure()
x0=4; y0=4;
width=3.5; height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
plot(2*pi./env.T, body.hydro.nu55{i,j}.*nu55_normalization,'-','LineWidth',1.5,...
    'Color','k','DisplayName','Analytical'); hold on
plot(Val.Data1.w,squeeze(Val.Data1.B(5,5,:))*Val.Data1.rho.*Val.Data1.w.',...
    '--','Color','k','DisplayName','WAMIT')
ylabel(nu55_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim(omega_limits)
%legend('Interpreter','Latex','Location','Best','NumColumns',2)
set(gca,'FontSize',12)

% combined mu55 and nu55
figure(); hold on; 
box on
x0=4; y0=4;
width=4.5; height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
% plot mu
yyaxis left
    % Create dummy legend entries:
    han = zeros(2,1);
    han(1) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.25,'DisplayName','$\mu_{55}$');hold on
    han(2) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor',nu_color,'MarkerEdgeColor',nu_color,'LineWidth',1.25,'DisplayName','$\nu_{55}$');
    %
    mu_color = 'k';
    plot(env.omega*omega_normalization, ...
        body.hydro.mu55{i,j}*mu55_normalization,'-','LineWidth',1.5,...
        'Color',mu_color,'DisplayName','Analytical');hold on
    plot(Val.Data1.w,squeeze(Val.Data1.A(5,5,:))*Val.Data1.rho,...
        '--','LineWidth',1.25,'Color',mu_color,'DisplayName','WAMIT')
    ylabel(mu55_label,'Interpreter','Latex');
    xlabel(omega_label,'Interpreter','Latex')
    xlim(omega_limits)
    set(gca,'FontSize',12)
% plot nu
yyaxis right
    nu_color = [135, 7, 7]/255;
    plot(2*pi./env.T, body.hydro.nu55{i,j}.*nu55_normalization,'-','LineWidth',1.5,...
        'Color',nu_color,'DisplayName','Analytical'); hold on
    plot(Val.Data1.w,squeeze(Val.Data1.B(5,5,:))*Val.Data1.rho.*Val.Data1.w.',...
        '--','LineWidth',1.25,'Color',nu_color,'DisplayName','WAMIT')
    ylabel(nu55_label,'Interpreter','Latex');
    xlabel(omega_label,'Interpreter','Latex')
    xlim(omega_limits)
    set(gca,'FontSize',12)
ax=gca; ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
legend(han,'Interpreter','Latex','Location','NorthEast')

% Added mass 15
figure()
x0=4; y0=4;
width=3.5; height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
plot(env.omega*omega_normalization, ...
    body.hydro.mu15{i,j}*mu15_normalization,'-','LineWidth',1.5,...
    'Color','k','DisplayName','Analytical');hold on
plot(Val.Data1.w,squeeze(Val.Data1.A(1,5,:))*Val.Data1.rho,...
    '--','LineWidth',1.5,'Color','k','DisplayName','WAMIT')
ylabel(mu15_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim(omega_limits)
%legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',12)

% radiation damping 15
figure()
x0=4; y0=4;
width=3.5; height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
plot(2*pi./env.T, body.hydro.nu15{i,j}.*nu15_normalization,'-','LineWidth',1.5,...
    'Color','k','DisplayName','Analytical'); hold on
plot(Val.Data1.w,squeeze(Val.Data1.B(1,5,:))*Val.Data1.rho.*Val.Data1.w.',...
    '--','LineWidth',1.5,'Color','k','DisplayName','WAMIT')
ylabel(nu15_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim(omega_limits)
ylim([0 850])
%legend('Interpreter','Latex','Location','Best','NumColumns',2)
set(gca,'FontSize',12)

% combined mu15 and nu15
figure(); hold on;
box on
x0=4; y0=4;
width=4.5; height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
% plot mu
yyaxis left
    % Create dummy legend entries:
    han = zeros(2,1);
    han(1) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.25,'DisplayName','$\mu_{15}$');hold on
    han(2) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor',nu_color,'MarkerEdgeColor',nu_color,'LineWidth',1.25,'DisplayName','$\nu_{15}$');
    %
    mu_color = 'k';
    plot(env.omega*omega_normalization, ...
        body.hydro.mu15{i,j}*mu15_normalization,'-','LineWidth',1.5,...
        'Color',mu_color,'DisplayName','Analytical');hold on
    plot(Val.Data1.w,squeeze(Val.Data1.A(1,5,:))*Val.Data1.rho,...
        '--','LineWidth',1.25,'Color',mu_color,'DisplayName','WAMIT')
    ylabel(mu15_label,'Interpreter','Latex');
    xlabel(omega_label,'Interpreter','Latex')
    xlim(omega_limits)
    set(gca,'FontSize',12)
% plot nu
yyaxis right
    nu_color = [135, 7, 7]/255;
    plot(2*pi./env.T, body.hydro.nu15{i,j}.*nu15_normalization,'-','LineWidth',1.5,...
        'Color',nu_color,'DisplayName','Analytical'); hold on
    plot(Val.Data1.w,squeeze(Val.Data1.B(1,5,:))*Val.Data1.rho.*Val.Data1.w.',...
        '--','LineWidth',1.25,'Color',nu_color,'DisplayName','WAMIT')
    ylabel(nu15_label,'Interpreter','Latex');
    xlabel(omega_label,'Interpreter','Latex')
    xlim(omega_limits)
    set(gca,'FontSize',12)
ax=gca; ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
legend(han,'Interpreter','Latex','Location','NorthEast')

% Linear restoring coefficient
figure
colororder({'k','#50F'})
plot(env.omega*omega_normalization, ...
    body.hydro.C55*ones(size(env.omega)),'LineWidth',1.25,...
    'DisplayName','Analytical Linear Restoring Coefficient');hold on
ylabel(C55_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
xlim(omega_limits)
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
xlim(omega_limits)
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
xlim(omega_limits)
legend('Interpreter','Latex','Location','SouthEast')
set(gca,'FontSize',16)


% Excitation torque magnitude and phase
Phase_color = [135, 7, 7]/255;
figure()
x0=4; y0=4;
width=4; height=3;
width=4.5; 
set(gcf,'units','inches','position',[x0,y0,width,height])
yyaxis left
    % Create dummy legend entries:
    han = zeros(2,1);
    han(1) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.25,'DisplayName','$|X_5|$');hold on
    han(2) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor',Phase_color,'MarkerEdgeColor',Phase_color,'LineWidth',1.25,'DisplayName','$\varphi_5$');
    
    plot(env.omega*omega_normalization,abs(body.hydro.X5{1,1})*X5_ma_normalization/1000,'-k',...
        'LineWidth',1.5,'DisplayName','Analytical Excitation Moment Magnitude'); hold on
    plot(Val.Data1.w,squeeze(Val.Data1.ex_ma(5,1,:))*Val.Data1.rho.*Val.Data1.g/1000,...
        '--','LineWidth',1.25,'Color','k','DisplayName',[Val.Source1,' Excitation Moment Magnitude'])
    xlabel(omega_label,'Interpreter','Latex')
    ylabel(X5_ma_label,'Interpreter','Latex')
    xlim(omega_limits);  ytickformat('%.1f')
    ylim([0 3.25])
    set(gca,'FontSize',12)
    legend(han,'Interpreter','Latex','Location','NorthEast')
yyaxis right
    plot(env.omega*omega_normalization,pi-angle(body.hydro.X5{1,1})*X5_ph_normalization,'-',...
        'Color',Phase_color,'LineWidth',1.5,'DisplayName','Analytical'); hold on
    plot(Val.Data1.w,squeeze(Val.Data1.ex_ph(5,1,:)),...
        '--','LineWidth',1.25,'Color',Phase_color ,'DisplayName',[Val.Source1])
    xlabel(omega_label,'Interpreter','Latex')
    ylabel(X5_ph_label,'Interpreter','Latex')
    ylim([0 1.6]); ytickformat('%.1f')
    set(gca,'FontSize',12)
ax=gca; ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
legend(han,'Interpreter','Latex','Location','NorthEast')


% Excitation force magnitude and phase
Phase_color = [135, 7, 7]/255;
figure()
x0=4; y0=4;
width=4; height=3.5;
width=4.5; height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
yyaxis left
    % Create dummy legend entries:
    han = zeros(2,1);
    han(1) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.25,'DisplayName','$|X_1|$');hold on
    han(2) = plot(nan,nan,'s','MarkerSize',5,'MarkerFaceColor',Phase_color,'MarkerEdgeColor',Phase_color,'LineWidth',1.25,'DisplayName','$\varphi_1$');
    plot(env.omega*omega_normalization,abs(body.hydro.X1{1,1})*X1_ma_normalization/1000,'-k',...
        'LineWidth',1.5,'DisplayName','Analytical Excitation Moment Magnitude'); hold on
    plot(Val.Data1.w,squeeze(Val.Data1.ex_ma(1,1,:))*Val.Data1.rho.*Val.Data1.g/1000,...
        '--','LineWidth',1.5,'Color','k','DisplayName',[Val.Source1,' Excitation Moment Magnitude'])
    xlabel(omega_label,'Interpreter','Latex')
    ylabel(X1_ma_label,'Interpreter','Latex')
    xlim(omega_limits);  ytickformat('%.1f')
    set(gca,'FontSize',12)
    legend(han,'Interpreter','Latex','Location','NorthEast')
yyaxis right
    plot(env.omega*omega_normalization,pi-angle(body.hydro.X1{1,1})*X1_ph_normalization,'-',...
        'Color',Phase_color,'LineWidth',1.5,'DisplayName','Analytical'); hold on
    plot(Val.Data1.w,squeeze(Val.Data1.ex_ph(1,1,:)),...
        '--','LineWidth',1.5,'Color',Phase_color ,'DisplayName',[Val.Source1])
    xlabel(omega_label,'Interpreter','Latex')
    ylabel(X1_ph_label,'Interpreter','Latex')
    ylim([0 1.6]); ytickformat('%.1f')
    set(gca,'FontSize',12)
ax=gca; ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
legend(han,'Interpreter','Latex','Location','NorthEast')



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



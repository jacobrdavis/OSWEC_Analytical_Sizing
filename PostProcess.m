%---------------------------- Header -------------------------------------%
% Title: VGOSWEC Sizing Tool General Post Processing Script
% Author: J. Davis
% Date: 12-15-20
% Last Update: 12-15-20
%-------------------------------------------------------------------------%
close all

% Visual results
normalize = false;          % Nondimensionalize plot results?
plot_indices = [6,1];       % Visualize output of individual case

%% SURFACE PLOTS
nrows = length(param.c_list);
ncols = length(param.w_list);

if ncols > 1 || nrows > 1 
%     figure
%     surf(body.dim.w,body.dim.c,out.maxCWR)
%     xlabel('Width [m]')
%     ylabel('Distance From Seabed [m]')
%     zlabel('CWR')
% 
%     figure
%     surf(body.dim.w,body.dim.c,out.maxFr1)
%     xlabel('Width [m]')
%     ylabel('Distance From Seabed [m]')
%     zlabel('Fr1')
    

%     figure
%     surf(body.dim.c,env.omega,[out.Xr1{:,1}])
%     xlabel('Distance From Seabed [m]')
%     ylabel('Angular Frequency \omega [rad/s]')
%     zlabel('Xr1')
%     
%     figure
%     surf(body.dim.c,env.omega,[out.CWR{:,1}])
%     xlabel('Distance From Seabed [m]')
%     ylabel('Angular Frequency \omega [rad/s]')
%     zlabel('CWR')

    if length(param.w_list) == 1 ||  length(param.c_list) == 1 
        xlab = 'Angular Frequency $\omega [rad/s]$';
        ylab = 'Distance From Seabed $[m]$';
        
        [plotout,axesout] = stackedlines(env.omega,body.dim.c,[out.Fr1{:,1}],true);
        zlab = 'Fr1';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
        
        if solver.calculateHydroEff == true  
            [plotout,axesout] = stackedlines(env.omega,body.dim.c,[out.HydroEff{:,1}]*100,true);
            zlab = 'Hydrodynamic Efficiency $[\%]$';
            xlabel(xlab,'Interpreter','Latex');
            ylabel(ylab,'Interpreter','Latex');
            zlabel(zlab,'Interpreter','Latex');           
        end
        
        [plotout,axesout] = stackedlines(env.omega,body.dim.c,[out.CWR{:,1}],true);
        zlab = 'CWR';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
    end
    
    if solver.calculateACE == 1
        figure
        surf(body.dim.w,body.dim.c,out.ACCW)
        xlabel('Width [m]')
        ylabel('Distance From Seabed [m]')
        zlabel('Avg ACCW')
        
        figure
        surf(body.dim.w,body.dim.c,cellfun(@(x) x(2),out.ACE))
        xlabel('Width [m]')
        ylabel('Distance From Seabed [m]')
        zlabel('ACE Med [m/$M]')
        
    end
    
end

%% INDIVIDUAL RESULTS PLOTS
if ncols==1 && nrows==1 
    i=1; j=1; 
else
    i = plot_indices(1); j = plot_indices(2);
end

body.dim.c(i,j)
% Normalization
if normalize == 1
    % Kelly 2017 normalizations:
%     mu_normalization = (rho*h^3*w*t)^(-1); mu_label = '\mu_{55}* = \mu_{55}/(\rhoh^3wt) [-]'; %[Kg-m^2]/[Kg-m^2]
%     nu_normalization = (rho*g^(1/2)*h^(5/2)*w*t)^(-1); nu_label = '\nu_{55}* = \nu_{55}/(\rhog^{1/2}h^{5/2}wt) [-]'; %[Kg-m^2/s]/[Kg-m^2/s]
    % Tom 2017 normalizations:
    mu_normalization = (body.prop.I55(i,j))^(-1); mu_label = '\mu_{55}* = \mu_{55}/(I_{55}) [-]'; %[Kg-m^2]/[Kg-m^2]
    nu_normalization = (env.omega.*body.prop.I55(i,j)).^(-1); nu_label = '\nu_{55}* = \nu_{55}/(\omegaI_{55}) [-]'; %[Kg-m^2/s]/[Kg-m^2/s]
    nu_g_normalization = (env.omega.*body.prop.I55(i,j)).^(-1); nu_g_label = '\nu_g* = \nu_g/(\omegaI_{55}) [-]'; %[Kg-m^2/s]/[Kg-m^2/s]
    Cg_normalization = (env.omega.^2.*body.prop.I55(i,j)).^(-1); Cg_label = 'C_g* = C_g/(\omega^2I_{55}) [-]'; %[Kg-m^2/s]/[Kg-m^2/s]
    omega_normalization = sqrt(env.h/env.g); omega_label = '\omega* = \omega (h/g)^{1/2} [-]'; % [rad/s]/[rad/s]
    X5_normalization = (env.rho*env.g*env.h^2*body.dim.w(i,j))^(-1); X5_label = '|X_5|* = |X5|/(\rhogh^2w) [-]'; %[N]/[N]
else
    mu_normalization = 1; mu_label = '\mu_{55} [Kg-m^2]';
    nu_normalization = 1; nu_label = '\nu_{55} [Kg-m^2/s]';
    nu_g_normalization = 1; nu_g_label = 'PTO Damping \nu_g [Kg-m^2/s]';
    Cg_normalization = 1; Cg_label = 'PTO Restoring Coefficient C_g [Kg-m^2/s^2]';
    omega_normalization = 1; omega_label = 'Angular Frequency \omega [rad/s]';
    X5_normalization = 1; X5_label = 'Excitation Torque |X_5| [N]';
end

% Added mass and radiation damping
figure
colororder({'k','#50F'})
yyaxis left;plot(env.omega*omega_normalization, ...
    body.hydro.mu55{i,j}*mu_normalization,'LineWidth',1.25);hold on
ylabel(mu_label);
xlabel(omega_label)
xlim([min(env.omega) max(env.omega)]*omega_normalization)
% ylimright = ylim;
yyaxis right;plot(2*pi./env.T, body.hydro.nu55{i,j}.*nu_normalization,'LineWidth',1.25);
ylabel(nu_label);
if exist('muVal','var') == 1 
plot(wVal,muVal.mu,'--','Color',colororder(1),...
    'DisplayName',[Val.Source,' Added Mass'])
end
if exist('nuVal','var') == 1 
plot(wVal,nuVal.nu,'--','Color','#50F',...
    'DisplayName',[Val.Source,' Added Mass'])
end

% ylim([ylimright])

% Excitation torque
figure
plot(env.omega*omega_normalization,body.hydro.X5{i,j}*X5_normalization,'k')
xlabel(omega_label)
ylabel(X5_label)
xlim([min(env.omega) max(env.omega)]*omega_normalization)



% % Pitch RAO
% figure
% plot(env.omega.*omega_normalization,RAO,'k')
% xlabel(omega_label)
% ylabel('RAO')
% xlim([min(omega) max(omega)]*omega_normalization)

% P_T max under motion constraint
% figure
% plot(omega,TAP)
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('P_T max, constrained')
% xlim([min(omega) max(omega)])

% Comparison of PTO damping coeff
% figure
% plot(omega,nu_g_constrained,'-'); hold on;
% plot(omega,nu_g,'--')
% plot(omega,nu_g2,'--')
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('PTO Damping \nu_g [Kg-m^2/s]')
% xlim([min(omega) max(omega)])
% legend('Optimized, constrained','Michele Output','\nu_{55}*Coefficient')

% PTO damping and restoring coefficients
figure
colororder({'#50F','#004C99'})
yyaxis left; plot(env.omega,body.pto.nu_g{i,j},'LineWidth',1.25); hold on
ylabel('\nu_g [Kg-m^2/s]')
if ~strcmp(body.pto.ctrltype,'free')==1
ylim([min([min(body.pto.nu_g{i,j}) min(body.pto.Cg{i,j})])...
    1.1*max([max(body.pto.nu_g{i,j}) max(body.pto.Cg{i,j})])])
end

yyaxis right; plot(env.omega,body.pto.Cg{i,j},'LineWidth',1.25); hold on
xlabel('Angular Frequency \omega [rad/s]')
ylabel('C_g [Kg-m^2/s^2]') %PTO Restoring Coefficient C_g [Kg-m^2/s^2]
xlim([min(env.omega) max(env.omega)])
if ~strcmp(body.pto.ctrltype,'free')==1
ylim([min([min(body.pto.nu_g{i,j}) min(body.pto.Cg{i,j})])...
    1.1*max([max(body.pto.nu_g{i,j}) max(body.pto.Cg{i,j})])])
end

% PTO spring-to-damping ratio
figure
plot(env.omega,body.pto.Cg{i,j}./(env.omega.*body.pto.nu_g{i,j}),'LineWidth',1.25)
ylabel('PTO Spring-to-Damping Ratio C_g/(\omega\nu_g) [-]')
xlabel('Angular Frequency \omega [rad/s]')
xlim([min(env.omega) max(env.omega)])

% Cf_opt output from Michele
% figure
% plot(omega,Cf_opt)
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('Cf opt from Michele')
% xlim([min(omega) max(omega)])

% delta
% figure
% plot(omega,delta)
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('\delta')
% xlim([min(omega) max(omega)])

% % TAP coeff and Vg/k
% figure
% yyaxis left
% plot(omega,TAP_Coeff); hold on
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('TAP Coefficient [-]')
% %xlim([min(omega) max(omega)])
% xlim([0.2 max(omega)])
% 
% yyaxis right
% plot(omega,Vg./k)
% ylabel('Vg/k [m^2/s]')
% ylim([0 300])

% % TAP/A^2
% figure
% plot(omega,TAPdivA2)
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('TAP per Wave Amplitude Squared P_T/A^2 [W/m^2]')
% xlim([min(omega) max(omega)])

% torque = (body.prop.I55(i,j) + body.hydro.mu55{i,j}).*env.omega.^2.*out.xi5{i,j}
% thetadot = 180/(2*pi)*env.omega.*out.xi5{i,j}


% Time-averaged PTO and wave power
figure
plot(env.omega,out.TAP{i,j}); hold on;
plot(env.omega,env.TAPwave*body.dim.w(i,j))
legend('PTO TAP P_T','Wave TAP w*P_w')
xlim([min(env.omega) max(env.omega)])
xlabel('Angular Frequency \omega [rad/s]')
ylabel('Power [W]')

% % Time-averaged PTO and wave power
% figure
% plot(env.T,out.TAP{i,j}); hold on;
% plot(env.T,env.TAPwave*body.dim.w(i,j))
% legend('PTO TAP P_T','Wave TAP w*P_w')
% xlim([min(env.T) max(env.T)])
% xlabel('Angular Frequency \omega [rad/s]')
% ylabel('Power [W]')

% Capture width ratio
figure
plot(env.omega,out.CWR{i,j})
xlim([min(env.omega) max(env.omega)])
xlabel('Angular Frequency $\omega [rad/s]$','Interpreter','Latex')
ylabel('Capture Width [m/m]')
set(gca,'FontSize',16)


% Surge force
F1_label = 'F1 [N]';
figure
plot(env.omega*omega_normalization,body.hydro.F1{i,j},'k')
xlabel(omega_label)
ylabel(F1_label)
xlim([min(env.omega) max(env.omega)]*omega_normalization)

% Hinge reaction forces
figure
plot(env.omega,out.Fr1{i,j}); hold on
xlim([min(env.omega) max(env.omega)])
xlabel('Angular Frequency $\omega [rad/s]$','Interpreter','Latex')
ylabel('Hinge Surge Reaction Force [N]')
set(gca,'FontSize',16)
if exist('Fr1Val','var') == 1 
plot(wVal,Fr1Val,'^',...
    'DisplayName',[Val.Source,' Fr1'])    
end

% Hydro efficiency
if solver.calculateHydroEff == true
    % Hydro Efficiency
    figure
    plot(env.omega,out.rotKE{i,j}); hold on
    plot(env.omega,out.projAreaWaveKEperW{i,j}*body.dim.w(i,j))
    xlim([min(env.omega) max(env.omega)])
    xlabel('Angular Frequency $\omega [rad/s]$','Interpreter','Latex')
    ylabel('Kinetic Energy [J]')
    set(gca,'FontSize',16)
end

% RAO
figure
plot(env.omega,out.RAO{i,j},'DisplayName','Analytical'); hold on
xlim([min(env.omega) max(env.omega)])
xlabel('Angular Frequency $\omega [rad/s]$','Interpreter','Latex')
% ylabel('RAO $i\omega\xi_{5}/A$','Interpreter','Latex')
ylabel('RAO $\xi_{5}/A$','Interpreter','Latex')
set(gca,'FontSize',16)
if exist('RAOVal','var') == 1 
plot(wVal,RAOVal,'^',...
    'DisplayName',[Val.Source,' RAO'])    
end
legend()

% Pitch amp xi5
figure
plot(env.omega,out.xi5{i,j}*180/pi,'DisplayName','Analytical'); hold on
xlim([min(env.omega) max(env.omega)])
xlabel('Angular Frequency $\omega [rad/s]$','Interpreter','Latex')
% ylabel('RAO $i\omega\xi_{5}/A$','Interpreter','Latex')
ylabel('Pitch Amplitude $\xi_{5} [deg]$','Interpreter','Latex')
set(gca,'FontSize',16)
if exist('xi5Val','var') == 1 
plot(wVal,xi5Val*180/pi,'^',...
    'DisplayName',[strrep(Val.Source,'_',' ')])    
end
legend()

% ACE Plots
if solver.calculateACE == 1
figure

% ACE seastates TAP
barlabels = cell(1,height(env.seastates));
for ss = 1:height(env.seastates)
barlabels{1,ss} = ['Tp = ',num2str(env.seastates.Tp(ss)),' s','\newline',...
                  'Hs = ',num2str(env.seastates.Hs(ss)),' m'];
end    

X = categorical(barlabels);
X = reordercats(X,barlabels);
bar(X,out.ACETAP{i,j}*10^(-3))
disp(barlabels{1,1})
ylabel('Time-averaged power [kW]')
    
end

%% POSTER FIGURES
post_fig = false;
if post_fig == 1
close all
mu_label = '$ A_{55} \left[kg-m^2\right]$';
nu_label = '$ B_{55} \left[\frac{kg-m^2}{s}\right]$';
nu_g_label = 'PTO Damping $ B_{PTO} \left[\frac{kg-m^2}{s}\right]$';
Cg_label = 'PTO Restoring Coefficient $ C_{PTO} \left[\frac{kg-m^2}{s^2}\right]$';
omega_label = 'Angular Frequency $\omega \left[\frac{rad}{s}\right]$';
X5_label = 'Excitation Torque $|X_5| \left[N\right]$';

mu_label = 'Added Mass $ A_{55} \left[kg.m^2\right]$';
nu_label = 'Radiation Damping $ B_{55} \left[\frac{kg.m^2}{s}\right]$';
nu_g_label = 'PTO Damping $ B_{PTO} \left[\frac{kg.m^2}{s}\right]$';
Cg_label = 'PTO Restoring $ C_{PTO} \left[\frac{kg.m^2}{s^2}\right]$';
omega_label = 'Angular Frequency $\omega \left[\frac{rad}{s}\right]$';
X5_label = 'Excitation Torque $|X_5| \left[\frac{N.m}{m}\right]$';


% Added mass and radiation damping
figure
colororder({'k','k'})
yyaxis left; plot(omega*omega_normalization, mu55*mu_normalization,'-k','LineWidth',2);hold on
ylabel(mu_label,'Interpreter','latex');
xlabel(omega_label,'Interpreter','latex')
xlim([min(omega) max(omega)]*omega_normalization)
% ylimright = ylim;

yyaxis right;plot(2*pi./T, nu55.*nu_normalization,'--k','LineWidth',2);
ylabel(nu_label,'Interpreter','latex');
set(gca,'FontSize',16)
legend('$A_{55}$','$B_{55}$','Interpreter','latex')
% ylim([ylimright])

% Excitation torque
figure
plot(omega*omega_normalization,X5*X5_normalization,'-k','LineWidth',2)
xlabel(omega_label,'Interpreter','latex')
ylabel(X5_label,'Interpreter','latex')
xlim([min(omega) max(omega)]*omega_normalization)
set(gca,'FontSize',16)




% PTO damping and restoring coefficients
figure
colororder({'k','k'})
yyaxis left; plot(omega,nu_g.*nu_g_normalization,'-k','LineWidth',2); hold on
ylabel(nu_g_label,'Interpreter','latex')

yyaxis right; plot(omega,Cg.*Cg_normalization,'--k','LineWidth',2); hold on
xlabel(omega_label,'Interpreter','latex')
ylabel(Cg_label,'Interpreter','latex') %PTO Restoring Coefficient C_g [Kg-m^2/s^2]
xlim([min(omega) max(omega)])
legend('$B_{PTO}$','$C_{PTO}$','Interpreter','latex')
set(gca,'FontSize',16)

% Time-averaged PTO and wave power
figure
plot(omega,TAP*10^-3,'-k','LineWidth',2); hold on;
plot(omega,TAPwave*w*10^-3,'--','Color',[0.0265 0.6137 0.8135],'LineWidth',2)
legend('PTO TAP $P_T$','Wave TAP $wP_w$','Interpreter','latex')
xlim([min(omega) max(omega)]); ylim([0 8])
xlabel(omega_label,'Interpreter','latex')
ylabel('Power $\left[kW\right]$','Interpreter','latex')
set(gca,'FontSize',16)

% Capture width ratio
figure
plot(omega,CWR,':r','LineWidth',2)
xlim([min(omega) max(omega)]); ylim([ 0 0.5])
xlabel(omega_label,'Interpreter','latex')
ylabel('CWR','Interpreter','latex')
set(gca,'FontSize',16)



end



%% SUB-FUNCTIONS

function [plotout, axesout] = stackedlines(x,y,zMat,transp_overlap)

if isrow(x); x = x.'; end

xMat = repmat(x, 1, length(y)); %// For plot3

%// Define y values
if ~isrow(y); y = y.'; end
yMat = repmat(y, numel(x), 1); %//For plot3

plotout = figure; axesout = gca
plot3(xMat, yMat, zMat, 'k','LineWidth',1); %// Make all traces blue
grid;
view(40,40); %// Adjust viewing angle so you can clearly see data
xlim([min(x) max(x)])
ylim([min(y) max(y)])
axis tight

if exist('transp_overlap', 'var') && transp_overlap == true
    ZL = zlim(gca);
    DZ = 0.03*(ZL(2)-ZL(1));
    
    for k=1:size(xMat,2)
        hPatch(k) = patch( ...
            [xMat(:,k);    flipud(xMat(:,k))   ], ...
            [yMat(:,k);    flipud(yMat(:,k))   ], ...
            [zMat(:,k);    flipud(zMat(:,k))-DZ], ...
            'w');
        set(hPatch(k), 'EdgeColor','none', 'FaceColor','w', 'FaceAlpha',0.9 );
    end
end

end

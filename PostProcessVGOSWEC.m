%---------------------------- Header -------------------------------------%
% Title: VGOSWEC Sizing Tool General Post Processing Script
% Author: J. Davis
% Date: 12-15-20
% Last Update: 12-15-20
%-------------------------------------------------------------------------%
%close all

% Visual results
normalize = false;          % Nondimensionalize plot results?
plot_indices = [2,1];       % Visualize output of individual case

% Define normalizations and labels:
omega_normalization = 1;    omega_label = '$\omega$ (rad s\textsuperscript{-1})';
% omega_normalization = sqrt(env.h/env.g); omega_label = '$\omega* = \omega (h/g)^{1/2} [-]$'; % [rad/s]/[rad/s]
mu55_normalization    = 1;	mu55_label    = '$\mu_{55}$ (kg-m\textsuperscript{2})';
mu15_normalization    = 1;	mu15_label    = '$\mu_{15}$ (kg-m)';
nu55_normalization    = 1;	nu55_label    = '$\nu_{55}$ (kg-m\textsuperscript{2}s\textsuperscript{-1})';
nu15_normalization    = 1;	nu15_label    = '$\nu_{15}$ (kg-m s\textsuperscript{-1})';
X5_ma_normalization = 1;	X5_ma_label = '$|X_5|$ (N-m)';
X5_ph_normalization = 1;	X5_ph_label = '$\varphi_5$ (rad)';
X1_ma_normalization = 1;	X1_ma_label = '$|X_1|$ (N)';
X1_ph_normalization = 1;	X1_ph_label = '$\varphi_1$ (rad)';
C55_normalization   = 1;	C55_label   = '$C_{55}$ (Kgm\textsuperscript{2}s\textsuperscript{-2})';
xi5_normalization   = 1;	xi5_label   = 'Pitch Amplitude $\xi_{5} [deg]$';
RAO_normalization   = 1;	RAO_label   = 'RAO $\xi_{5}/A$';
F1_ma_normalization   = 1;	F1_ma_label    = '$|F_1|$ (N)';
F5_ma_normalization = 1;	F5_ma_label    = '$|F_5|$ (N-m)';
Fr1_normalization   = 1;	Fr1_label   = 'Hinge Surge Reaction Force Magnitude $[N]$';
Cg_normalization   = 1;     Cg_label   = '$C_{g}$ (Kgm\textsuperscript{2}s\textsuperscript{-2})';
nu_g_normalization    = 1;	nu_g_label    = '$\nu_{g}$ (kg-m\textsuperscript{2}s\textsuperscript{-1})';
Lnorm = (env.h)^(-1);

%     figure
%         x = 0.01;
%         y = 0.0001;
%         loglog(x,y,'x')
%         xtickvals = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2] ; xticks(xtickvals); 
%         ytickvals = [0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05]; yticks(ytickvals); 
%         xlims = [0.001 0.2]; xlim(xlims)
%         ylims = [0.00005 0.05] ; ylim(ylims)
%         
%         
%         pbaspect([1 695/596 1])
%         axis tight;
%         hold on
%         I = imread('Water_wave_theories.png');
%         h = image([0 1],[ 0 1],I);
%         uistack(h,'bottom')
%         % https://www.mathworks.com/matlabcentral/answers/96023-how-do-i-add-a-background-image-to-my-gui-or-figure-window

%figure
% hold on
% plot(repmat(env.T,ncols,nrows),abs([body.hydro.F1{:,1}]))
% 
% zlab = F1_ma_label;
% xlabel('T (s)','Interpreter','Latex');
% ylabel(zlab,'Interpreter','Latex');
% steepness_list(i,1) = env.steepness;
% i = i+1
% legend([repmat('h/L = ',length(steepness_list),1),num2str(steepness_list)])

%% SURFACE PLOTS
nrows = length(param.c_list);
ncols = length(param.w_list);

% % availible power ma
% figure
% plot(repmat(env.omega,ncols,nrows),0.5*abs([body.hydro.F5{:,1}]).*(env.omega.*abs([out.xi5{:,1}])))
% legend([[repmat('I_{55} = ',nrows,ncols),num2str(body.prop.I55,2),repmat(' kg-m^2',nrows,ncols)],....
%         [repmat('; C_{55} = ',nrows,ncols),num2str(body.hydro.C55,2),repmat(' kg-m^2 s^-2',nrows,ncols)]])
% zlab = 'Available Power (W)';
% xlabel('Angular Frequency $\omega$ (rad/s)','Interpreter','Latex');
% ylabel(zlab,'Interpreter','Latex');
% 
% figure
% plot(repmat(env.T,ncols,nrows),.5*abs([body.hydro.F5{:,1}]).*(env.omega.*abs([out.xi5{:,1}])))
% legend([[repmat('I_{55} = ',nrows,ncols),num2str(body.prop.I55,2),repmat(' kg-m^2',nrows,ncols)],....
%         [repmat('; C_{55} = ',nrows,ncols),num2str(body.hydro.C55,2),repmat(' kg-m^2 s^-2',nrows,ncols)]])
% zlab = 'Available Power (W)';
% xlabel('T (s)','Interpreter','Latex');
% ylabel(zlab,'Interpreter','Latex');

% % availible power ma
% figure
% plot(repmat(env.omega,ncols,nrows),0.5*abs([body.hydro.F5{:,1}]).*(env.omega.*abs([out.xi5{:,1}])))
% legend([[repmat('t = ',nrows,ncols),num2str(body.dim.t,2),repmat(' m',nrows,ncols)],....
%         [repmat('; C_{55} = ',nrows,ncols),num2str(body.hydro.C55,2),repmat(' kg-m^2 s^-2',nrows,ncols)]])
% zlab = 'Available Power (W)';
% xlabel('Angular Frequency $\omega$ (rad/s)','Interpreter','Latex');
% ylabel(zlab,'Interpreter','Latex');
% 
% figure
% plot(repmat(env.T,ncols,nrows),.5*abs([body.hydro.F5{:,1}]).*(env.omega.*abs([out.xi5{:,1}])))
% legend([[repmat('t = ',nrows,ncols),num2str(body.dim.t,2),repmat(' m',nrows,ncols)],....
%         [repmat('; C_{55} = ',nrows,ncols),num2str(body.hydro.C55,2),repmat(' kg-m^2 s^-2',nrows,ncols)]])
% zlab = 'Available Power (W)';
% xlabel('T (s)','Interpreter','Latex');
% ylabel(zlab,'Interpreter','Latex');

if ncols > 1 || nrows > 1 

%-------------------------------------------------------------------------%
% frequency (\omega) and distance to seabed (c) plots
%-------------------------------------------------------------------------%
    if length(param.w_list) == 1 
        xlab = 'Angular Frequency $\omega [rad/s]$';
        ylab = 'Distance From Seabed $[m]$';
        
        % added mass
        figure
        plot(repmat(env.omega,ncols,nrows),[body.hydro.mu55{:,1}])
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = mu55_label;
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');

        figure
        plot(repmat(env.T,ncols,nrows),[body.hydro.mu55{:,1}])
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = mu55_label;
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        % rad damping
        figure
        plot(repmat(env.omega,ncols,nrows),[body.hydro.nu55{:,1}])
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = nu55_label;
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
        plot(repmat(env.T,ncols,nrows),[body.hydro.mu55{:,1}])
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = nu55_label;
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        % ex ma 
        figure
        plot(repmat(env.omega,ncols,nrows),abs([body.hydro.X5{:,1}]))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = X5_ma_label;
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
        plot(repmat(env.T,ncols,nrows),abs([body.hydro.X5{:,1}]))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = X5_ma_label;
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
         % torque ma 
        figure
        plot(repmat(env.omega,ncols,nrows),abs([body.hydro.F5{:,1}]))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = F5_ma_label;
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
        plot(repmat(env.T,ncols,nrows),abs([body.hydro.F5{:,1}]))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = F5_ma_label;
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
         % availible power ma 
        figure
        plot(repmat(env.omega,ncols,nrows),0.5*abs([body.hydro.F5{:,1}]).*(env.omega.*abs([out.xi5{:,1}])))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = 'Available Power (W)';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
        plot(repmat(env.T,ncols,nrows),.5*abs([body.hydro.F5{:,1}]).*(env.omega.*abs([out.xi5{:,1}])))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = 'Available Power (W)';
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        % Wave amp
        figure
        plot(env.omega,env.Aw)
        zlab = 'Wave Amp (m)';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
        plot(env.T,env.Aw)
        zlab = 'Wave Amp (m)';
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        text(2,0.005,['H/L = ',num2str(env.steepness)])
        2*env.Aw(5)/(env.g*env.T(5)^2)
        2*env.Aw(1)/env.lambda(1)
        
        % pitch amp
        figure
        plot(repmat(env.omega,ncols,nrows),180/pi*abs([out.xi5{:,1}]))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = xi5_label;
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
        plot(repmat(env.T,ncols,nrows),180/pi*abs([out.xi5{:,1}]))
        legend([repmat('c = ',nrows,ncols),num2str(body.dim.c),repmat(' m',nrows,ncols)])
        zlab = xi5_label;
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        
        [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.c,abs([out.Fr1{:,1}]),true);
        zlab = 'Fr1';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
        
        if solver.calculateHydroEff == true  
            [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.c,[out.HydroEff{:,1}]*100,true);
            zlab = 'Hydrodynamic Efficiency $[\%]$';
            xlabel(xlab,'Interpreter','Latex');
            ylabel(ylab,'Interpreter','Latex');
            zlabel(zlab,'Interpreter','Latex');           
        end
        
        [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.c,[out.CWR{:,1}],true);
        zlab = 'CWR';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
        
        if solver.calculateACE == 1
            figure
            plot(body.dim.c,out.ACCW,'x')
            xlabel('Distance From Seabed $[m]$','Interpreter','Latex');
            ylabel('Avg ACCW [m]','Interpreter','Latex');

            figure
            plot(body.dim.c,cellfun(@(x) x(2),out.ACE),'x')
            xlabel('Distance From Seabed $[m]$','Interpreter','Latex');
            ylabel('ACE Med [m/$M]','Interpreter','Latex');
            
            figure
            plot(body.dim.c,cellfun(@(x) x(2),out.CCE),'x')
            xlabel('Distance From Seabed $[m]$','Interpreter','Latex');
            ylabel('CCE Med [$M]','Interpreter','Latex');
            
        end
        
%-------------------------------------------------------------------------%
% frequency (\omega) and width (w) plots
%-------------------------------------------------------------------------%
    elseif length(param.c_list) == 1 
        xlab = 'Angular Frequency $\omega$ [rad/s]';
        ylab = 'Width [m]';
       
         % availible power ma 
        figure
        plot(repmat(env.omega,nrows,ncols),0.5*abs([body.hydro.F5{1,:}]).*(env.omega.*abs([out.xi5{1,:}])))
        legend([repmat('w = ',nrows,ncols),num2str(body.dim.w),repmat(' m',nrows,ncols)])
        zlab = 'Available Power (W)';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        figure
         plot(repmat(env.T,nrows,ncols),0.5*abs([body.hydro.F5{1,:}]).*(env.omega.*abs([out.xi5{1,:}])))
        legend([repmat('w = ',ncols,nrows),num2str(body.dim.w.'),repmat(' m',ncols,nrows)])
        zlab = 'Available Power (W)';
        xlabel('T (s)','Interpreter','Latex');
        ylabel(zlab,'Interpreter','Latex');
        
        %
        
        figure
        plot(repmat(env.omega,nrows,ncols),[body.hydro.mu55{1,:}])
        legend([repmat('w = ',ncols,nrows),num2str(transpose(body.dim.w)),repmat(' m',ncols,nrows)])
        
        
        [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.w,abs([body.hydro.mu55{1,:}]),true);
        zlab = '$\mu_{55}$ [kg-m\textsuperscript{2}]';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
        
        [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.w,abs([out.Fr1{1,:}]),true);
        zlab = '$F_{r1}$ [N]';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
        
        if env.defspectrum == false
        [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.w,abs([out.TAP{1,:}]),true);
        zlab = '$P_T$ [W]';
        xlabel(xlab,'Interpreter','Latex');
        ylabel(ylab,'Interpreter','Latex');
        zlabel(zlab,'Interpreter','Latex');
        end
        
        if solver.calculateHydroEff == true  
            [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.w,[out.HydroEff{1,:}]*100,true);
            zlab = 'Hydrodynamic Efficiency [%]';
            xlabel(xlab,'Interpreter','Latex');
            ylabel(ylab,'Interpreter','Latex');
            zlabel(zlab,'Interpreter','Latex');           
        end
        
        if env.defspectrum == false
            [plotout,axesout] = fun.graphics.stackedlines(env.omega,body.dim.w,[out.CWR{1,:}],true);
            zlab = 'CWR [m/m]';
            xlabel(xlab,'Interpreter','Latex');
            ylabel(ylab,'Interpreter','Latex');
            zlabel(zlab,'Interpreter','Latex');
        end
        
        if solver.calculateACE == 1
            figure
            plot(body.dim.w,out.ACCW,'x')
            xlabel('Width [m]','Interpreter','Latex');
            ylabel('Avg ACCW [m]','Interpreter','Latex');

            figure
            plot(body.dim.w,cellfun(@(x) x(2),out.ACE),'x')
            xlabel('Width [m]','Interpreter','Latex');
            ylabel('ACE Med [m/\$M]','Interpreter','Latex');
            
            figure
            plot(body.dim.w,cellfun(@(x) x(2),out.CCE),'x')
            xlabel('Width [m]','Interpreter','Latex');
            ylabel('CCE Med [\$M]','Interpreter','Latex');
            
        end    
%-------------------------------------------------------------------------%
% width (w) and distance to seabed (c) plots 
%-------------------------------------------------------------------------%
    elseif length(param.w_list) > 1 &&  length(param.c_list) > 1
 
        xlab = '$w/h$';
        ylab = '$c/h$';
        xlimits = [round(min(body.dim.w(1,:)*Lnorm),2) max(body.dim.w(1,:)*Lnorm)];
        ylimits = [min(body.dim.c(:,1)*Lnorm) max(body.dim.c(:,1)*Lnorm)];
        xtickvals =0.:0.1:0.9;
        ytickvals =0:0.1:.6;
        surfcolor = 'jet'; % bone; jet
        fontsize = 12;
        
        if solver.calculateACE == 1
            
            % ACCW mcp plot
            figure; hold on
            % zlimits = [0.55 0.7]; zlim(zlimits) 
            ztickform = '%g'; zlab = 'ACCW (m)';
            mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm ,out.ACCW);
            run('mcpsettings.m')
            
            % CCE mcp plot
            figure; hold on
            % zlimits = [0.55 0.7]; zlim(zlimits) 
            ztickform = '%g'; zlab = 'CCE Med (\$M)';
            mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) x(2),out.CCE));
            run('mcpsettings.m')
            
            % ACE mcp plot
            figure; hold on
            % zlimits = [0.55 0.7]; zlim(zlimits) 
            ztickform = '%g'; zlab = 'ACE Med (m/\$M)';
            mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) x(2),out.ACE));
            run('mcpsettings.m')
              
        end
         
      
        % Capture width ratio mcp plot       
        figure; hold on
        %zlimits = [0.55 0.7]; zlim(zlimits); 
        ztickform = '%.2f'; zlab = 'CWR (-)';
        mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) max(x),out.CWR));
        run('mcpsettings.m')

        % TAP mcp plot       
        figure; hold on
        %zlimits = [0.55 0.7]; zlim(zlimits); 
        ztickform = '%.2f'; zlab = '$P_{WEC}$ (W)';
        mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) max(x),out.TAP));
        run('mcpsettings.m')
        
        % Fdn D mcp plot
        figure; hold on
        % zlimits = [0.55 0.7]; zlim(zlimits)
        ztickform = '%g'; zlab = 'Foundation Diameter (m)';
        mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,fdn.dim.D);
        run('mcpsettings.m')
        
        % Fdn D mcp plot
        figure; hold on
        % zlimits = [0.55 0.7]; zlim(zlimits)
        ztickform = '%g'; zlab = '$F_{R1}$ (MN)';
        mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,abs(out.maxFr1)*10^(-6));
        run('mcpsettings.m')
        
        % Fdn Mb mcp plot
        figure; hold on
        % zlimits = [0.55 0.7]; zlim(zlimits)
        ztickform = '%g'; zlab = '$M_{b}$ (N-m)';
        mcp = meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,out.MbFdnBaseMax);
        run('mcpsettings.m')
        
    end
end

%% INDIVIDUAL RESULTS PLOTS
indv = 0
if indv == true
if ncols==1 && nrows==1 
    i=1; j=1; 
else
    i = plot_indices(1); j = plot_indices(2);
end



% Pitch Added mass and radiation damping
figure
yyaxis left
plot(env.omega*omega_normalization, ...
    body.hydro.mu55{i,j}*mu55_normalization,'k','LineWidth',1.25,...
    'DisplayName','Added Mass'); hold on
ylabel(mu55_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex');
xlim([min(env.omega) max(env.omega)]*omega_normalization)
% ylimright = ylim;
yyaxis right
plot(2*pi./env.T, body.hydro.nu55{i,j}.*nu55_normalization,'k:','LineWidth',1.25,...
    'DisplayName','Radiation Damping');
ylabel(nu55_label,'Interpreter','Latex');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',16)

% Surge-Pitch Added mass and radiation damping
figure
yyaxis left
plot(env.omega*omega_normalization, ...
    body.hydro.mu15{i,j}*mu15_normalization,'k','LineWidth',1.25,...
    'DisplayName','Added Mass'); hold on
ylabel(mu15_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex');
xlim([min(env.omega) max(env.omega)]*omega_normalization)
% ylimright = ylim;
yyaxis right
plot(2*pi./env.T, body.hydro.nu15{i,j}.*nu15_normalization,'k:','LineWidth',1.25,...
    'DisplayName','Radiation Damping');
ylabel(nu15_label,'Interpreter','Latex');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',16)

% Linear restoring coefficient
figure
plot(env.omega*omega_normalization, ...
    body.hydro.C55(i,j)*ones(size(env.omega))*C55_normalization,'k','LineWidth',1.25);hold on
xlim([min(env.omega) max(env.omega)]*omega_normalization)
ylabel(C55_label,'Interpreter','Latex');
xlabel(omega_label,'Interpreter','Latex')
set(gca,'FontSize',16)

% Excitation torque magnitude and phase
figure
yyaxis left
plot(env.omega*omega_normalization,abs(body.hydro.X5{i,j})*X5_ma_normalization,'k',...
    'DisplayName','Magnitude'); hold on
xlabel(omega_label,'Interpreter','Latex')
ylabel(X5_ma_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)

yyaxis right
plot(env.omega*omega_normalization,pi-angle(body.hydro.X5{i,j})*X5_ph_normalization,'k:',...
    'DisplayName','Phase'); hold on
xlabel(omega_label,'Interpreter','Latex')
ylabel(X5_ph_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',16)


% Excitation force magnitude and phase
figure
yyaxis left
plot(env.omega*omega_normalization,abs(body.hydro.X1{i,j})*X1_ma_normalization,'k',...
    'DisplayName','Magnitude'); hold on
xlabel(omega_label,'Interpreter','Latex')
ylabel(X1_ma_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)

yyaxis right
plot(env.omega*omega_normalization,pi-angle(body.hydro.X1{i,j})*X1_ph_normalization,'k:',...
    'DisplayName','Phase'); hold on
xlabel(omega_label,'Interpreter','Latex')
ylabel(X1_ph_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',16)

% PTO damping and restoring coefficients
figure
yyaxis left; plot(env.omega,body.pto.nu_g{i,j}*nu_g_normalization,'k','LineWidth',1.25,...
    'DisplayName','Damping'); hold on
ylabel(nu_g_label,'Interpreter','Latex')
if ~strcmp(body.pto.ctrltype,'free')==1
ylim([min([min(body.pto.nu_g{i,j}) min(body.pto.Cg{i,j})])...
    1.1*max([max(body.pto.nu_g{i,j}) max(body.pto.Cg{i,j})])])
end

yyaxis right; plot(env.omega,body.pto.Cg{i,j}*Cg_normalization,'k:','LineWidth',1.25,...
    'DisplayName','Restoring'); hold on
xlabel(omega_label,'Interpreter','Latex')
ylabel(Cg_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)])
if ~strcmp(body.pto.ctrltype,'free')==1
ylim([min([min(body.pto.nu_g{i,j}) min(body.pto.Cg{i,j})])...
    1.1*max([max(body.pto.nu_g{i,j}) max(body.pto.Cg{i,j})])])
end
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
legend('Interpreter','Latex','Location','Best')
set(gca,'FontSize',16)



% Time-averaged PTO and wave power
figure
plot(env.omega,out.TAP{i,j}); hold on;
plot(env.omega,env.TAPwave*body.dim.w(i,j))
legend('PTO TAP P_T','Wave TAP w*P_w')
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
ylabel('Power [W]')
set(gca,'FontSize',16)

% Capture width ratio
figure
plot(env.omega,out.CWR{i,j})
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
ylabel('Capture Width [m/m]')
set(gca,'FontSize',16)


% Surge force
figure
plot(env.omega*omega_normalization,abs(body.hydro.F1{i,j}),'k')
xlabel(omega_label,'Interpreter','Latex')
ylabel(F1_ma_label,'Interpreter','Latex')
xlim([min(env.omega) max(env.omega)]*omega_normalization)
set(gca,'FontSize',16)

% Hinge reaction forces
figure
plot(env.omega,abs(out.Fr1{i,j}),'DisplayName','Surge'); hold on
plot(env.omega,abs(out.Fr3{i,j}),'DisplayName','Heave');
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
ylabel('Hinge Reaction Force [N]')
set(gca,'FontSize',16)


% Hydro efficiency
if solver.calculateHydroEff == true
    % Hydro Efficiency
    figure
    plot(env.omega,out.rotKE{i,j}); hold on
    plot(env.omega,out.projAreaWaveKEperSA{i,j}*body.dim.w(i,j))
    xlim([min(env.omega) max(env.omega)])
    xlabel(omega_label,'Interpreter','Latex')
    ylabel('Kinetic Energy [J]')
    set(gca,'FontSize',16)
end

% RAO
figure
plot(env.omega,abs(out.RAO{i,j}),'DisplayName','Analytical'); hold on
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
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
plot(env.omega,abs(out.xi5{i,j}*180/pi),'DisplayName','Analytical'); hold on
xlim([min(env.omega) max(env.omega)])
xlabel(omega_label,'Interpreter','Latex')
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


body.dim.w(i,j); body.dim.t(i,j); body.dim.ht(i,j)

Do = 2.58;

gamma = 0.75

Di = Do*gamma

Cc = Do/2

Icy = pi/64*(Do^4-Di^4)
Mb = max(abs(out.Fr1{i,j}))*body.dim.c(i,j)
SigB= Mb*Cc/Icy % Pa
Sigyield = 70E6 % Pa
  
Sigyield/SigB

SF = 2
Do = (Mb/(pi/32*Sigyield/SF*(1-gamma^4)))^(1/3)
end

%% SUB-FUNCTIONS

% function [plotout, axesout] = stackedlines(x,y,zMat,transp_overlap)
% 
% if isrow(x); x = x.'; end
% 
% xMat = repmat(x, 1, length(y)); %// For plot3
% 
% %// Define y values
% if ~isrow(y); y = y.'; end
% yMat = repmat(y, numel(x), 1); %//For plot3
% 
% plotout = figure; axesout = gca
% plot3(xMat, yMat, zMat, 'k','LineWidth',1); %// Make all traces blue
% grid;
% view(40,40); %// Adjust viewing angle so you can clearly see data
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
% axis tight
% 
% if exist('transp_overlap', 'var') && transp_overlap == true
%     ZL = zlim(gca);
%     DZ = 0.03*(ZL(2)-ZL(1));
%     
%     for k=1:size(xMat,2)
%         hPatch(k) = patch( ...
%             [xMat(:,k);    flipud(xMat(:,k))   ], ...
%             [yMat(:,k);    flipud(yMat(:,k))   ], ...
%             [zMat(:,k);    flipud(zMat(:,k))-DZ], ...
%             'w');
%         set(hPatch(k), 'EdgeColor','none', 'FaceColor','w', 'FaceAlpha',0.9 );
%     end
% end
% 
% end

%         % mesh setttings
%         mcp(1).EdgeColor = 'flat'; 
%         mcp(1).FaceColor = 'flat';
%         mcp(1).FaceAlpha = 0.7;
%         mcp(1).FaceLighting ='gouraud';
%         % contour settings
%         mcp(2).EdgeColor = 'flat'; 
%         xlabel(xlab,'Interpreter','Latex','FontSize',fontsize,'Position',[0.6785,-0.0550,0.5340])
%         ylabel(ylab,'Interpreter','Latex','FontSize',fontsize,'Position',[0.2515, 0.3371,0.5361])
%         zlabel('CWR [-]','Interpreter','Latex','FontSize',fontsize)
%         xlim(xlimits)
%         ylim(ylimits)
%         xticks(xtickvals); 
%         ztickformat(ztickform)
%         yticks(ytickvals)
%         view(-30,30)
%         ax = gca
%         ax.XTickLabelRotation=-35.5;
%         ax.YTickLabelRotation=12;

%         figure
%         colormap(surfcolor)
%         surf(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) max(x),out.TAP),out.maxCWR,...
%             'FaceColor','flat','FaceAlpha',0.75,...
%             'EdgeColor','flat','FaceLighting','gouraud')
%         xlabel(xlab,'Interpreter','Latex')
%         ylabel(ylab,'Interpreter','Latex')
%         zlabel('TAP [W]','Interpreter','Latex')
%         c = colorbar;
%         c.Label.String = 'CWR [m/m]';
%         xlim(xlimits)
%         ylim(ylimits)
%         xticks(xtickvals)
%         yticks(ytickvals)
%         
%         figure
%         colormap(surfcolor)
%         surf(body.dim.w*Lnorm,body.dim.c*Lnorm,out.maxCWR,cellfun(@(x) max(x),out.TAP),...
%             'FaceColor','flat','FaceAlpha',0.75,...
%             'EdgeColor','flat','FaceLighting','gouraud')
%         xlabel(xlab,'Interpreter','Latex')
%         ylabel(ylab,'Interpreter','Latex')
%         zlabel('CWR [-]','Interpreter','Latex')
%         c = colorbar
%         c.Label.String = 'TAP [W]'
%         xlim(xlimits)
%         ylim(ylimits)
%         xticks(xtickvals)
%         yticks(ytickvals)
     
%         figure; hold on
%         x0=4; y0=4;
%         width=4.5; height=3.5;
%         set(gcf,'units','inches','position',[x0,y0,width,height])
%         colormap(surfcolor)
%         surf(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) max(x),out.CWR),...
%             'FaceColor','flat','FaceAlpha',0.7,...
%             'EdgeColor','flat','FaceLighting','gouraud')
%        figure
%         meshc(body.dim.w*Lnorm,body.dim.c*Lnorm,cellfun(@(x) max(x),out.CWR))
%         
%         planeimg = abs(cellfun(@(x) max(x),out.CWR))
%         zLimits = get(gca,'ZLim')
%         imgzposition =zLimits(1); % desired z position of the image plane.
%        
%         surf(xlimits,ylimits,repmat(imgzposition, [2 2]),...
%             planeimg,'facecolor','texture')
% 
%         % set the view angle.
% %         view(45,30);
%         
%         
%         xlabel(xlab,'Interpreter','Latex','FontSize',fontsize)
%         ylabel(ylab,'Interpreter','Latex','FontSize',fontsize)
%         zlabel('CWR [-]','Interpreter','Latex','FontSize',fontsize)
%         xlim(xlimits)
%         ylim(ylimits)
%         xticks(xtickvals)
%         yticks(ytickvals)  

%
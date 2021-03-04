% Uniform mesh countour plot settings script for EWTEC '21 plots

% size and position
x0=4; y0=4;
width=4.5; height=3.5;
width=3.5; height=3.5*3.5/4.5;
% xlabpos = [0.6785,-0.0550,0.5340]; ylabpos = [0.2515, 0.3371,0.5361];
ylabpos = [0.106665545222148,0.083208933438449 - 0.03,0];
xlabpos = [0.725983912702808,0.02154291981955 - 0.03,0];
set(gcf,'units','inches','position',[x0,y0,width,height])
% axes graphics
box on; grid on;
% color
colormap(surfcolor)
% mesh setttings
mcp(1).EdgeColor = 'flat';
mcp(1).FaceColor = 'flat';
mcp(1).FaceAlpha = 0.7;
mcp(1).FaceLighting ='gouraud';
% contour settings
mcp(2).EdgeColor = 'flat';
% labels, limits, and ticks
xlabel(xlab,'Interpreter','Latex','FontSize',fontsize,'units','norm','Position',xlabpos)
ylabel(ylab,'Interpreter','Latex','FontSize',fontsize,'units','norm','Position',ylabpos)
zlabel(zlab,'Interpreter','Latex','FontSize',fontsize)
xlim(xlimits)
ylim(ylimits)
xticks(xtickvals);
ztickformat(ztickform)
yticks(ytickvals)
% view and label orientation
view(-30,30)
ax = gca;
ax.XTickLabelRotation=-35.5;
ax.YTickLabelRotation=12;
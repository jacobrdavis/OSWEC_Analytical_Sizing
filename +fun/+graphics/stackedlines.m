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

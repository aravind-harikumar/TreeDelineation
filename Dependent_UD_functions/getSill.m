function sill = getSill(csmNormalized, isPlotOn)
if(isPlotOn)
    figure;
    subplot(2,2,1);
    imshow(csmNormalized);
end
[cent, varargout] = FastPeakFind(mat2gray(csmNormalized),10);
getX = cent(1:2:end);
getY = cent(2:2:end);

x = getX;
y = getY;
z = sub2ind(size(csmNormalized),getY,getX);

if(isPlotOn)
    hold on;
    scatter(x,y,4,z,'filled'); box on;
    ylabel('y'); xlabel('x')
    %h = labelpoints (x, y, num2cell(z), 'Color', 'w'); % labels points
    title('data (coloring according to z-value)')
    
    subplot(2,2,2);
    hist(z,20);
    ylabel('frequency'); xlabel('z');
    title('histogram of z-values');

    subplot(2,2,3);
end

d = variogram([x y],z,'plotit',isPlotOn,'nrbins',50);

if(isPlotOn)
    hold on;
end

d.val  = ((d.val-min(d.val))/(min(d.val)-max(d.val)));
d.val = abs(d.val);
c0 = max(d.val); a0 = max(d.distance)*2/3; %a0 = 15; c0 = 0.1;
[range,sill,nugget] = variogramfit(d.distance,d.val,a0,c0,[],'model','exponential', 'solver','fminsearchbnd', 'nugget',0, 'plotit',isPlotOn);

if(isPlotOn)
    title('Isotropic variogram');
    subplot(2,2,4);
end

d2 = variogram([x y],z,'plotit',isPlotOn,'nrbins',50,'anisotropy',true);
if(isPlotOn)
    title('Anisotropic variogram');
end
end
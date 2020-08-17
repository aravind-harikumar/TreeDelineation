function plotPlaneFrom3Points(pointA,pointB,pointC, isMirrored)

    normal = cross(pointA-pointB, pointA-pointC); %# Calculate plane normal
    %# Transform points to x,y,z
    x = [pointA(1) pointB(1) pointC(1)];  
    y = [pointA(2) pointB(2) pointC(2)];
    z = [pointA(3) pointB(3) pointC(3)];

    %Find all coefficients of plane equation    
    A = normal(1); B = normal(2); C = normal(3);
    D = -dot(normal,pointA);
    %Decide on a suitable showing range
    xLim = [min(x)*-10 max(x)*10];
    zLim = [min(z)*-10 max(z)*10];
    [X,Z] = meshgrid(xLim,zLim);
    Y = (A * X + C * Z + D)/ (-B);
    reOrder = [1 2  4 3];
    hold on;
    if(isMirrored)
        patch(X(reOrder),Y(reOrder),Z(reOrder),'k');
        patch(-X(reOrder),-Y(reOrder),Z(reOrder),'k');
    else
        patch(X(reOrder),Y(reOrder),Z(reOrder),'k');
    end
    grid on;
    alpha(0.3);
    
end
% 
% function plotPlaneFrom3Points(P1,P2,P3, plotXScale, plotYScale, C)
%     normal = cross(P1-P2, P1-P3);
%     if(normal(3) == 0)
%         normal(3) = 1; 
%     end
% 
%     Xn= normal(1); Yn= normal(2); Zn= normal(3);
%     Px= P1(1); Py= P1(2); Pz= P1(3);     
%     [X, Y] = meshgrid(plotXScale, plotYScale);
%     
%     Fn = @(Xn,Yn,Zn,Px,Py,Pz,X,Y)(-Xn*X + Xn*Px - Yn*Y + Yn*Py + Zn*Pz)/Zn;
%     Z = Fn(Xn,Yn,Zn,Px,Py,Pz,X,Y);
% 
%     h(2)=surf(X, Y, Z);
%     %set(h,'edgecolor','k');
%     set(h(2), 'edgecolor','none')
%     set(h(2),'CData',C);
%     %alpha(h(2),0.5)
%     %colormap(summer);
%     
%     %syms x y z
%     %P = [x,y,z];
%     %planefunction = dot(normal, P-P1);
%     %zplane = solve(planefunction, z);
%     %fmesh(zplane,  [-1, 50, -1, 10]);
% end
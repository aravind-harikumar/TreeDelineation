function plotShape(coeffArr, xmid, ymid, maxZVal, isParabloid, pltParaboloid)
    a = coeffArr(1); b = coeffArr(2); c = coeffArr(3);
    
    X1= [-2.5:0.1:2.5];
    Y1= [-2.5:0.1:2.5];
    [X,Y] = meshgrid(X1,Y1);
    
    Z = [];
    if(isParabloid)
        % Plot Parabloid
        Z = myparabloid(X,Y,a,b,c);  
    else
        % Plot myCone
        Z = myCone(X,Y,a,b,c);  
    end 
    if(pltParaboloid)
        %subplot(1,2,1);
        surf(X+xmid,Y+ymid,Z+maxZVal);
    end
end
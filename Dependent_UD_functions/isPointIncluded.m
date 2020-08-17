function isIncluded = isPointIncluded(point,a,b,c,r,shapeType)   
    isIncluded = 0;%a=1.4;b=1.4;c=12.4;
    x = point(1); y = point(2); z = point(3);
    Z =[];
    if (strcmp(shapeType, 'elcone'));
        Z = myCone(x,y,a,b,c);
    elseif (strcmp(shapeType, 'elparaboloid'));
        Z = myparabloid(x,y,a,b,c);
    elseif(strcmp(shapeType,'elCylinder'));        
        dis = sqrt(sum((x - y).^ 2));
        if(dis > r)
            Z = z+1;
        else
            Z = z-1;
        end
    end
    if(z>Z)
        isIncluded = 1;
    elseif(z==Z)
        isIncluded = 2;
    end
end
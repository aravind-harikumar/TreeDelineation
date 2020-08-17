function getDist = getDist2Point(point1, point2, extPoint)
    x = extPoint; %some point
    a = point1; %segment points a,b
    b = point2;
    
    v = b-a;
    w = x-a;
    c1 = w*v';
    c2 = v*v';
    if(c1<=0)
       getDist =  sqrt(sum((x - a).^ 2)); 
       return;
    elseif (c2<=c1)
       getDist =  sqrt(sum((x - b).^ 2)); 
       return;
    else
        b1= c1/c2;
        pb= a + b1*v;
        getDist =  sqrt(sum((x - pb).^ 2));   
    end
  
end
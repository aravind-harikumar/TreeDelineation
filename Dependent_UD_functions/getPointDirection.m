function signVal = getPointDirection(ppa, ppb, ppc, externPoint) 
   [a,b,c,d] = getPlaneEqParameters(ppa, ppb, ppc);   
    signVal = sign(a*externPoint(1) + b*externPoint(2) + c*externPoint(3) + d);
end

function [a,b,c,d] = getPlaneEqParameters(ppa, ppb, ppc)  
    ppab = ppb - ppa;
    ppac = ppc - ppa;
    nVec = cross(ppab,ppac);   
    a = nVec(1); b = nVec(2); c = nVec(3);
    d = -(a*ppa(1) + b*ppa(2) + c*ppa(3)); 
end
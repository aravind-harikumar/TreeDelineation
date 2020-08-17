function avgVal = getRC(r,c,a)
    if(and(and(r>1,r<size(a,1)),and(c>1,c<size(a,2))))
        %avgVal = (a(r,c)+a(r+1,c)+a(r-1,c)+a(r,c-1)+a(r,c+1))/4;
        avgVal = (a(r,c)+a(r+1,c)+a(r-1,c)+a(r,c-1)+a(r,c+1)+a(r+1,c+1)+a(r-1,c-1)+a(r+1,c-1)+a(r-1,c+1))/8;
    else
        avgVal = 0;
    end
end
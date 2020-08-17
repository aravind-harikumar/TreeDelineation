 function resultant = getVectorResultant(X,Y)
    %X(1) = X(1)+13; X(2) = X(2)-123;
    UNITX = X/norm(X); UNITY = Y/norm(Y);
    resultant = (UNITX - UNITY)/2;
 end
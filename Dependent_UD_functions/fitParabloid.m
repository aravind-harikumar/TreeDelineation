function coeffVal  = fitParabloid(inputPoints)
    x =inputPoints(:,1); y = inputPoints(:,2); z = inputPoints(:,3);

    x = x- mean(x);
    y = y- mean(y);

    A = [x.*x y.*y];
    
    isParabloid = true;    
    if(isParabloid)
        coeff =  A\z; % for hyperboloid
        a = coeff(1); b = coeff(2);
        c = mean(z(:)./ ( (power(x(:),2)*a^2) + (power(y(:),2)*b^2) ) ); %old 
        
%         if(a>b)
%             temp=b;
%             a=b;
%             a=temp;
%         end
%         if(c<-3)
%             c = -3;
%         end
%         %a = 21; b = 10; c = 0.05;
%         %a = ev(1); b = ev(2); c = 60;
%         %c=0.2;
% 
%         if(abs(a/b) > 2)
%         a= 25; b=10; c =-0.05;
%         elseif(abs(b/a) > 2)
%         a= 25; b=10; c =-0.05;
%         end
         a= 25; b=25; c =-0.04;
        
    else 
        coeff =  A\(z.*z);  % for cone
        a = coeff(1); b = coeff(2);
        c = mean(sqrt( power(z(:),2)./ ( (power(x(:),2)*a^2) + (power(y(:),2)*b^2) ) )); % cone
        %a = 25; b = 10; c = 0.2;
    end

    coeffVal = [a b c];
end
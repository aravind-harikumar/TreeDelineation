function Norm=normalize(A)
    Norm_A=[];
    for i=1:size(A,2)
        Norm_col=(A(:,i)-min(A(:,i)))/(max(A(:,i))-min(A(:,i)));
        Norm_A=[Norm_A,Norm_col];
    end
    Norm=Norm_A;
end
function retEigenData = getBestEigenDirection(inputData)
    [a,~,c]=pca(inputData);
    pc = zeros(3,2);
    
    for i=1:size(a,2)
        pc(:,i) = a(:,i);%'*-c(i);
    end
    retEigPoints{1} = pc;
    retEigPoints{2} = c;    
    retEigenData = retEigPoints;    
end
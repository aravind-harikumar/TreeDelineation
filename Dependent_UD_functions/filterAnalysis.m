function filterAnalysis(csm)
    %close all;
    figure;
    %load('csm.mat');
    xImg = size(csm,1); yImg = size(csm,2);
    
    f = makeLMfilters();
    [Fille,Filbr,Filhe] = size(f); 
    arrFilterResponse = zeros(xImg+Fille-1,yImg+Filbr-1,Filhe);
    for i=1:36
            subplot(6,6,i);
            arrFilterResponse(:,:,i) = conv2(csm, f(:,:,i));
            imagesc(arrFilterResponse(:,:,i))
    end        
    
end
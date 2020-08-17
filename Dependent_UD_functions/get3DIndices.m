function outPutSubMatices = get3DIndices(data, grixDim)
     
    outPutSubMatices={};
    getSubMatricesbyIndex = {};
    %subplot(1,2,1);
    searchGrid = [3 3 3];
    
    psz = (2*(floor(searchGrid(1)/2)-1));
    if(psz == 0)
        psz = 1;
    end

    dataArrPadded = padarray(data,[psz psz psz],-1); % padded with -1;
    
    [X,Y,Z] = meshgrid(linspace(1,grixDim(1),grixDim(1)+(psz)*2),linspace(1,grixDim(2),grixDim(2)+(psz)*2),linspace(1,grixDim(3),grixDim(3)+(psz)*2));
    searchIndex = size(X,1)*size(X,2)*size(X,3);
    
    %axis([-(grixDim(1)+2) grixDim(1)+2 -(grixDim(2)+2) grixDim(2)+2 -(grixDim(3)+2) grixDim(3)+2]); axis equal;
   
    sss =[];
    idx = [];
    for i = 1:1:searchIndex  
        
        if(i==1518)
           sdf = 0; 
        end
        
        [Iadj, ~, ~] = neighbourND(i, [grixDim(1)+(psz)*2 grixDim(2)+(psz)*2 grixDim(3)+(psz)*2], [2 2 2]);
        Iadj = [Iadj i];
        
        sss = [sss length(Iadj)];        
        tempdata = dataArrPadded;
        
        
        thrsh =  (searchGrid(1)^2)*2 + ((searchGrid(1)-1)^2)*(searchGrid(1)-2) + 1 ;
        if( and(length(Iadj)==searchGrid(1)*searchGrid(2)*searchGrid(3), sum(tempdata(Iadj)==-1)<thrsh))  % 
            getSubMatricesbyIndex{i} = Iadj; % valid indice 
        else
            getSubMatricesbyIndex{i} = Iadj;
            idx = [idx i];  % invalid indice
        end 
    end    
    
    %max(sss)
    
    validIndices = setdiff([1:1:searchIndex],idx);
    for j = 1:1:size(validIndices,2)
 
%         subplot(1,2,1);
%         hold off;
%         tempdata2 = dataArrPadded;
%         tempdata2(getSubMatricesbyIndex{validIndices(j)})=-100;
%         scatter3(X(:), Y(:), Z(:),2000, tempdata2(:) ,'.','MarkerFaceColor',[0 0 0.5]); 
%         text(X(:)+0.1, Y(:)+0.1, Z(:)+0.1, num2str(dataArrPadded(:)));
%         axis square;
%         
%         subplot(1,2,2);
        indi = getSubMatricesbyIndex{validIndices(j)};
        ngVals = dataArrPadded(indi); 
        outPutSubMatices{j} = reshape(ngVals,[3 3 3]);

%         scatter3(X(indi), Y(indi), Z(indi),100,'MarkerFaceColor',[0 0 0.5]); 
%         text(X(indi)+0.1, Y(indi)+0.1, Z(indi)+0.1, num2str(dataArrPadded(indi)'));
%         axis equal;
%         axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:)) max(Z(:))]); 
%         pause(0.1);
    end 
    
end
function getIDsOfTrees()
    % read the data of all trees
    %[~,~,raw] = xlsread('C:\My_Files\1_PhD_Research\2_PhD_Publications\8_TGRS_Paper_3\PLOTS\ALL_TREES.xlsx');
    [~,~,raw] = xlsread('C:\Users\harikumar\Desktop\Excel Files\ALL_PLOTS_COMBINED_1.xlsx');
    for idx = 2:numel(raw)
        if isnumeric(raw{idx})
            raw{idx} = num2str(raw{idx});
        end
    end
    raw(cellfun('isempty',raw))={0};
    
    XY = str2double(raw(2:end,6:7));
    
    IDList1 = [];
    [Idx,D] = rangesearch(XY(1:196,:),XY(1:196,:),5);
    for i = 1:1:size(Idx,1)        
        temp = cell2mat(Idx(i));
        if(size(temp,2)>1)
            IDList1 = [IDList1 temp];
        end
    end
    IDList1 = unique(IDList1)';    
    
    IDList2 = [];
    [Idx,D] = rangesearch(XY(1:196,:),XY(1:196,:),inf);
    for i = 1:1:size(Idx,1)        
        temp = cell2mat(Idx(i));
        if(size(temp,2)>1)
            IDList2 = [IDList2 temp];
        end
    end
    IDList2 = unique(IDList2)';    
    IDList3 = setdiff(IDList2,IDList1)';
    
    % < 2.5
    PropCV = str2double(raw(IDList1,40:42));
    MeanV =  mean(PropCV,1); % 68.7158    5.8006   -0.6998
    
    % between 2.5 and 5
    SoACV = str2double(raw(IDList3,47:49));
    MeanVd = mean(SoACV(2:25,:),1); % 57.1477    5.6973    0.2090
    
    %  > 5
    %78.3496    6.7280   -1.0584
    
    
    d =0;
    %out = pdist2(XY,XY,'euclidean');
    
    %out(out==0) = inf;    
    %[ooo,idx]= min(out,[],2);
    
%     out = tril(out);
%     
%     out(out==0) = inf;
%     out(isnan(out)) = inf;
%     out(or(out<2.5,out>5)) = 0;
%     %out(out>5) = 0;
%     %out(out~=0) = inf;
%     
%     AllInx = [];
%     for ii =1:1:196
%        tmpCol = out(:,ii);
%        indxx = find(tmpCol ~= 0);
%        AllInx = [AllInx; indxx];
%     end
%     
%     indxReal = unique(AllInx);
%     
%     %[indx, bal] = find(out == 0);
    
end
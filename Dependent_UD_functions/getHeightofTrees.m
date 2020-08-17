function getHeightofTrees()
    % read the data of all trees
    %[~,~,raw] = xlsread('C:\My_Files\1_PhD_Research\2_PhD_Publications\8_TGRS_Paper_3\PLOTS\ALL_TREES.xlsx');
    [~,~,raw] = xlsread('C:\Users\harikumar\Desktop\Excel Files\plots1_PropCR.csv');
    for idx = 2:numel(raw)
        if isnumeric(raw{idx})
            raw{idx} = num2str(raw{idx});
        end
    end
    raw(cellfun('isempty',raw))={0};
    
    %get specieswise data    
    %uniqueSpecies = unique(raw(:,5));    
    uniqueSpecies = {'ar','la','ab','pc'};
    
    %for i = 1:length(uniqueSpecies) 
    
    cond1 = ~isnan(str2double(raw(:,7)));
    cond2 = ~isnan(str2double(raw(:,8)));        
    matData11 = str2double(raw(and(cond1,cond2), 7:8));
    
    cond11 = ~isnan(str2double(raw(:,17)));
    cond22 = ~isnan(str2double(raw(:,18)));  
    matData22 = str2double(raw(and(cond11,cond22), 17:18));
    
    opArr =[];
    
    for i = 1:1:size(matData11,1)        
        distances = pdist2(matData11(i,:), matData22);
        [minDistance, indexOfMinDistance] = min(distances);
         opArr = [opArr; str2double(raw(indexOfMinDistance, 20))];
    end
  
end

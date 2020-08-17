function getheightfromDBHBySpecies()
    % read the data of all trees
    %[~,~,raw] = xlsread('C:\My_Files\1_PhD_Research\2_PhD_Publications\8_TGRS_Paper_3\PLOTS\ALL_TREES.xlsx');
    [~,~,raw] = xlsread('C:\Users\harikumar\Desktop\Excel Files\TREES_ORIGINAL_DATA.xlsx');
    for idx = 1:numel(raw)
        if isnumeric(raw{idx})
            raw{idx} = num2str(raw{idx});
        end
    end
    raw(cellfun('isempty',raw))={0};
    
    %get specieswise data    
    %uniqueSpecies = unique(raw(:,5));    
    uniqueSpecies = {'ar','la','ab','pc'};
    
    for i = 1:length(uniqueSpecies)        
        isSpeciesRow = strcmp(raw(:,5), uniqueSpecies(i));
        speciesdata = raw(isSpeciesRow,:);
        
        %get species rows with height and DBH
        
        cond1 = ~isnan( str2double(speciesdata(:,10)));
        cond2 = ~isnan( str2double(speciesdata(:,6)) );        
        speciesdatawithheightanddbh = speciesdata(and(cond1,cond2), :);
        matData = [str2double(speciesdatawithheightanddbh(:,6)) str2double(speciesdatawithheightanddbh(:,10))]; % dbh and height
        matData = table(log(matData(:,1)),matData(:,2),'VariableNames',{'DBH','HEIGHT'});

        % get regression model
        lm = fitlm(matData,'HEIGHT~DBH')
        
        % ar = (y intercept)-25.292 +  14.13;   
        % la = (y intercept)-34.267 +  15.494  
        % ab = (y intercept)-32.941 +  15.772
        %
        %plotResiduals(lm)
        %Res = table2array(lm.Residuals);
        %h = archtest(Res(:,1));
    end    
    
end

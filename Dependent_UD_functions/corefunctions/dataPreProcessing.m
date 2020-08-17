function dataAtt = dataPreProcessing(singleTreeLiDARdata)
% Write normalized data to a table for performance improvement 
        dataAtt.origData = write2table(singleTreeLiDARdata);
       
        % make point density uniform across height.
        %dataAtt.origData = increasePointDensityToUDValue(dataAtt.origData,false); % isploton = false/true;
        
        buffer = (max(dataAtt.origData(:,3)) - min(dataAtt.origData(:,3)))*0.15;        
        dataAtt.origData = dataAtt.origData( dataAtt.origData(:,3)>(min(dataAtt.origData(:,3)+buffer)),:); 
        
        lidarDataArr = normalizeLiDARData(dataAtt.origData);
        %lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:);
        
        dataAtt.treeWidth =  max(lidarDataArr(:,1)); treeHeight=max(lidarDataArr(:,2)); % to calculate max wisth/breadth to set plot width
        dataAtt.maxTreeWidth = max(dataAtt.treeWidth, treeHeight); % Keep it as 5 if issues arise
        dataAtt.mintreeHeight = min(lidarDataArr(:,3));
        dataAtt.maxtreeHeight = max(lidarDataArr(:,3));        
        dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;        
        dataAtt.htDeduction = 0; %dataAtt.treeHeight*0.05; % to get rid of ground noise points
        
        % Crown height
        lidarDataArr = lidarDataArr(find(and(lidarDataArr(:,3) >= min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) <= max(lidarDataArr(:,3)))),:);
        dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
        dataAtt.maxCrownHeight = max(lidarDataArr(:,3));
        dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;
                           
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);
        
        % Identify the center point of the tree (in top view).
        [dataAtt.retMaxXYZ, dataAtt.retMaxXYZIndex] = findMaxHeightXY(lidarDataArr);
        
         dataAtt.origDataMaxXY = dataAtt.origData(dataAtt.retMaxXYZIndex,1:2);
        %opOffseted = getPCAVectorsFotTreeTop(lidarDataArr(:,1:3));
        
        
        %lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); 
        %lidarDataArr = lidarDataArr(find(lidarDataDensityArr(:,1) >= 0.0000),:);
        
        %lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); %get point in very close neightbourhood                
        index = 1:1:size(lidarDataArr,1); 
         
        % % Density priuned LiDAR data
        %selIndx = find(lidarDataDensityArr(:,1) >= 0.0000);        
        dataAtt.lidarDataArray =  lidarDataArr; %lidarDataArr(selIndx,:);                
        dataAtt.index = index; %index(selIndx);
        dataAtt.lidarDataDensityArr = zeros(size(lidarDataArr,1),1); %lidarDataDensityArr(selIndx,:);        
        dataAtt.lidarDataArray = [dataAtt.lidarDataArray dataAtt.index' dataAtt.lidarDataDensityArr];
        
        dataAtt.lidarDataArray(:,1) = dataAtt.lidarDataArray(:,1) - dataAtt.retMaxXYZ(1);
        dataAtt.lidarDataArray(:,2) = dataAtt.lidarDataArray(:,2) - dataAtt.retMaxXYZ(2);
        
        
        [dataAtt.retMaxXYZNew, dataAtt.retMaxPointIndexNew] = findMaxHeightXYNew(dataAtt.lidarDataArray);
        
        dataAtt.XMax = max(dataAtt.lidarDataArray(:,1));
        dataAtt.XMin = min(dataAtt.lidarDataArray(:,1));

        dataAtt.YMax = max(dataAtt.lidarDataArray(:,2));
        dataAtt.YMin = min(dataAtt.lidarDataArray(:,2));

        dataAtt.ZMax = max(dataAtt.lidarDataArray(:,3));
        dataAtt.ZMin = min(dataAtt.lidarDataArray(:,3));
        
       
        
       % dataAtt.lidarDataArray = dataAtt.lidarDataArray(~and((abs(dataAtt.lidarDataArray(:,2)).^2 +abs(dataAtt.lidarDataArray(:,1)).^2)<0.1, dataAtt.lidarDataArray(:,3)<=0.9*max(dataAtt.lidarDataArray(:,3))),:);   

end




function newLiDARArr = increasePointDensityToUDValue(lidarDataArr,isploton)
    numDivisions = 10;    
    % find layer with max density
    layerPointDensity = zeros(numDivisions,1); 
    mintreeHeight = min(lidarDataArr(:,3));
    maxtreeHeight = max(lidarDataArr(:,3));
    %[mintreeHeight,maxtreeHeight] = deal(num2cell([min(lidarDataArr(:,3)),max(lidarDataArr(:,3))]));    
    divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
    for i = 1:1:size(divsArr,2)-1
        layerPointDensity(i) = size(lidarDataArr(and(lidarDataArr(:,3) >= divsArr(i), lidarDataArr(:,3) < divsArr(i+1)),:),1);
    end
    [maxDensity,maxDensityIndex] = max(layerPointDensity);
    maxDensity = 2500;
    newLiDARArr = [];
    % increase density of remaing layers to the maxDensity
    for i = size(divsArr,2)-1:-1:1
        currentCloud = lidarDataArr(and(lidarDataArr(:,3) >= divsArr(i), lidarDataArr(:,3) < divsArr(i+1)),:);
        %if(~isequal(i,maxDensityIndex))
           %lackingPointCount = maxDensity - size(currentCloud,1);         
           %randindices = int32(round(rand(maxDensity*((size(divsArr,2))-i),1)*(size(currentCloud,1)-1),0))+1;
           randindices = int32(round(rand(maxDensity,1)*(size(currentCloud,1)-1),0))+1;
           newLiDARArr = [newLiDARArr; currentCloud(randindices,:)];
           %plot3(newLiDARArr(:,1),newLiDARArr(:,2), newLiDARArr(:,3),'.','MarkerSize', 10,'Color',[0 0.5 0]);
           %axis equal;
        %else
        %   newLiDARArr = [newLiDARArr; currentCloud];
        %end
    end
end

% function newLiDARArr = increasePointDensityToUDValue(lidarDataArr,isploton)
%     numDivisions = 10;
%     
%     % find layer with max density
%     layerPointDensity = zeros(numDivisions,1); 
%     mintreeHeight = min(lidarDataArr(:,3));
%     maxtreeHeight = max(lidarDataArr(:,3));
%     %[mintreeHeight,maxtreeHeight] = deal(num2cell([min(lidarDataArr(:,3)),max(lidarDataArr(:,3))]));    
%     divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
%     for i = 1:1:size(divsArr,2)-1
%         layerPointDensity(i) = size(lidarDataArr(and(lidarDataArr(:,3) >= divsArr(i), lidarDataArr(:,3) < divsArr(i+1)),:),1);
%     end
%     [maxDensity,maxDensityIndex] = max(layerPointDensity);
%     
%     newLiDARArr = []; 
%     % increase density of remaing layers to the maxDensity
%     for i = 1:1:size(divsArr,2)-1
%         currentCloud = lidarDataArr(and(lidarDataArr(:,3) >= divsArr(i), lidarDataArr(:,3) < divsArr(i+1)),:);
%         
%         if(~isequal(i,maxDensityIndex)) 
%            lackingPointCount = maxDensity - size(currentCloud,1); %+ (50*(size(divsArr,2)-i));
%            
%            %[~, percCover] = getPointDensityCover(currentCloud,i,isploton);
%            %lackingPointCount = int32(lackingPointCount + lackingPointCount*(1-percCover));
%            
%            randindices = int32(round(rand(lackingPointCount,1)*(size(currentCloud,1)-1),0))+1;
%            newLiDARArr = [newLiDARArr; currentCloud; currentCloud(randindices,:)];
%         else
%            newLiDARArr = [newLiDARArr; currentCloud];
%         end
%     end
%     
% end

% function [csm, percCover] = getPointDensityCover(currentCloud,icnt,isploton)
%     MMS = getMaxMins(currentCloud);
%     gXStep = abs(MMS.maxX - MMS.minX)/50;
%     gYStep = abs(MMS.maxY - MMS.minY)/50;
%     halfXGridSize = gXStep/2;       
%     halfYGridSize = gXStep/2;
%     xDimGrid = size(MMS.minX+halfXGridSize : gXStep : MMS.maxX-halfXGridSize,2);
%     yDimGrid = size(MMS.minY+halfYGridSize : gYStep : MMS.maxY-halfYGridSize,2);  
%     
%     csm = zeros(xDimGrid,yDimGrid);
%     % Generating csm and cdm
%     
%     cnt1 = 1;   ijArr = [];  ijALL = []; 
%     for i = MMS.minX+halfXGridSize:gXStep:MMS.maxX-halfXGridSize
%         cnt2 = 1;
%         for j = MMS.minY+halfYGridSize:gYStep:MMS.maxY-halfYGridSize               
%             cond1 = and(currentCloud(:,1)>i-halfXGridSize, currentCloud(:,1)<=i+halfXGridSize);
%             cond2 = and(currentCloud(:,2)>j-halfYGridSize, currentCloud(:,2)<=j+halfYGridSize);
%             indGridPoints = find(and(cond1,cond2));
%             ijALL = [ijALL; cnt1 cnt2]; 
%             if(size(indGridPoints,1)>0)  
%                 csm(cnt1,cnt2) = numel(currentCloud(indGridPoints,3));
%                 ijArr = [ijArr; cnt1 cnt2]; 
%             end
%             cnt2 = cnt2+1;
%         end
%         cnt1 = cnt1+1;
%     end
%     
%     k = boundary(ijArr(:,1),ijArr(:,2),0);
%     boundaryArr = [ijArr(k,1) ijArr(k,2)];
%     
%     [in,on] = inpolygon(ijALL(:,2), ijALL(:,1), ijArr(k,1), ijArr(k,2));
%     AllPoinInOn = or(in,on);
%     csm = ((csm-min(csm(:)))/(max(csm(:))-min(csm(:))))*255;
%     percCover = size(ijArr,1)/sum(AllPoinInOn);
%     
%     if(isploton)
%         subplot(3,3,icnt);
%         imagesc(csm); hold on;
%         plot(boundaryArr(:,2),boundaryArr(:,1),'-*', 'Color', 'k');
%         plot(ijALL(AllPoinInOn,1),ijALL(AllPoinInOn,2),'.', 'MarkerSize', 5, 'Color', 'r');
%         plot(ijArr(:,2),ijArr(:,1),'.', 'MarkerSize', 10, 'Color', 'g');
%         title(percCover);
%     end
%    
% end



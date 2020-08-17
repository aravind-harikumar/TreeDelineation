function [projLiDARdata, DiffDist, EV_ratio_thresh, EVals] = getProjectedLiDARdata(stData, scaleFactor)
    
    plotProjectedData = true;
    %plotXScale = 0:1:stData.maxtreeHeight+0.5;
    %plotYScale = -10:0.1:10;
    projLiDARdata = zeros(size(stData.lidarDataArray,1),8);
    
    maxXYZ = stData.retMaxXYZ;
    %retMaxPointIndex = stData.retMaxPointIndexNew;
    PCAplotOn = true;
    [~, ~, REFResultpoint, DiffDist, EV_ratio_thresh, EVals]  = getPCAVectors(stData,PCAplotOn); 
    %REFpoint(1) = REFpoint(1) - 3;
    opOffseted = [REFResultpoint(1:2,:) [10;10] [0;0] [0;0] [stData.lidarDataArray(size(stData.lidarDataArray,1),6)+1;stData.lidarDataArray(size(stData.lidarDataArray,1),6)+1] ];
    
    orginalDataXYZ = [stData.lidarDataArray(:,1:3); opOffseted(1,1:3)];
    refPointXYZ = repmat(opOffseted(1,1:3), size(orginalDataXYZ,1),1);
    dataArrLidAR = [stData.lidarDataArray(:,1:6);opOffseted(1,1:6)];

    [~, indexExtreme]= min(abs(sqrt(sum((orginalDataXYZ.^2-refPointXYZ.^2), 2))));
    %retVal = [dataArrLidAR(indexExtreme,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexExtreme,3)];
    
    indexTreeTop = stData.retMaxXYZIndex;
    %retTreeTop = [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexTreeTop,3)];
    
    % Get the plain points for performing calcualtons
    tempPPa = [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexTreeTop,2) 0];
    
    %%% another method to find reference point
    %     [optimalRefPoint,~] = getOptimalPlaneRefPoint(stData);
    %     if(~isempty(optimalRefPoint))
    %         %EV_ratio_thresh = median(valleyPointMean);
    %         if(optimalRefPoint(1)<0)
    %             optimalRefPoint = -optimalRefPoint;
    %         end
    %     end
    
    if(isempty(REFResultpoint))
        tempPPb = [dataArrLidAR(indexExtreme,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexExtreme,3)];  
    else
        tempPPb = [REFResultpoint(1,1) REFResultpoint(1,2) dataArrLidAR(indexExtreme,3)];
    end 
    tempPPc = [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexTreeTop,2) dataArrLidAR(indexTreeTop,3)];
    
    % Perform Projection
    sz = size(dataArrLidAR,1);
    for i = 1:1:sz          
        inputPoint = dataArrLidAR(i,1:3);            
        inputPointIndex = dataArrLidAR(i,6);  
        externPoint = inputPoint;        
        position = getPointDirection(tempPPa, tempPPb, tempPPc, externPoint);
        [sectorLength, disttocentre] = getSectorLength(tempPPb, inputPoint, [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexTreeTop,2)], [maxXYZ(1) maxXYZ(2)], scaleFactor); %% disttoCentre = Y, sectorlength = X,  inputPointZ = Z  
        projLiDARdata(i,1:8) = [inputPoint(3) position*sectorLength  disttocentre 0 inputPoint inputPointIndex];    
    end    
    projLiDARdata = real(projLiDARdata);

    % Get the points for making the 3D plane
    ppa = [0 0 projLiDARdata(indexExtreme,3)];
    ppb = [projLiDARdata(indexTreeTop,1) 0 projLiDARdata(indexTreeTop,3)];
    ppc = [projLiDARdata(indexExtreme,1) 0 projLiDARdata(indexExtreme,3)];    
    plainPointsAll = [ppa;ppb;ppc;ppa];
    
    % plot result
    if(and(NC.ISPLOTON,plotProjectedData))        
        f2 = figure;  
        set(f2, 'Position', [10 60 600 420])           
        hold on; 
        % plot projected lidar data
        % uu = plot3(projLiDARdata(:,1), projLiDARdata(:,2), projLiDARdata(:,3),'.');
        %uu = scatter(projLiDARdata(:,1), projLiDARdata(:,2),300, projLiDARdata(:,3),'.');
       % projLiDARdata(:,3) = projLiDARdata(:,3) +  
        uu = scatter3(projLiDARdata(:,1), projLiDARdata(:,2), projLiDARdata(:,3),300, projLiDARdata(:,3),'.');
        xlabel('Tree Height'); ylabel('Distance to the refernce cyclinder surface'); zlabel('Distance to Stem');
        %camproj perspective; rotate3d on; axis vis3d; 
        axis equal; axis on; grid on;  box on; view(-130,60);
        %set(gca, 'Zdir', 'reverse')
        axis([min(projLiDARdata(:,1)) max(projLiDARdata(:,1)) min(projLiDARdata(:,2)) max(projLiDARdata(:,2)) min(projLiDARdata(:,3)) max(projLiDARdata(:,3))]);
        daspect([1 0.1 1]);
        title('Projected LiDAR Data');
        % generate colormap
        T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
        x = [0 85 110 255];
        map = interp1(x/255,T,linspace(0,1,100));
        colormap(map);
        %C = uu.CData;        
        % plot points used for generating the plane
        hold on;
        plot3(plainPointsAll(:,1), plainPointsAll(:,2) ,plainPointsAll(:,3), '.', 'MarkerSize', 40,'Color', [1 0 0]);
        % plot plane  
        plotPlaneFrom3Points(ppa, ppb, ppc, false);
        text(ppa(1),ppa(2),ppa(3),'Projection of tree top on XY plane','HorizontalAlignment','right');
        text(ppc(1),ppc(2),ppc(3),'reference point','HorizontalAlignment','right');
        text(ppb(1),ppb(2),ppb(3),'Tree top', 'HorizontalAlignment','right');
    end
    projLiDARdata(:,4) = [1:1:size(stData.lidarDataArray,1)+1]';
    
    % remove reference point added
    projLiDARdata = projLiDARdata(1:end-1,:);
    
   % projLiDARdata = reducePointHeight(projLiDARdata);
    
    % genereating basic info of projected data
    projLiDARdata = projDataPreProcessing(projLiDARdata);
    
    
    
end



function projDataDensified = increasePointDensityToUDValuesss(PJData)
   PJData = projDataPreProcessing(PJData); 
   VDim = getVoxelDimension(PJData, 0.5);
   projDataDensified = getVoxelDensityAndNebrCellIndicesNewsss(PJData, VDim.rX, VDim.rY, VDim.rZ);
end

function projDataDensified = getVoxelDensityAndNebrCellIndicesNewsss(PCData, xedge, yedge, zedge)

xlen = length(xedge)-1; xmid = zeros(xlen,1);
ylen = length(yedge)-1; ymid = zeros(ylen,1);
zlen = length(zedge)-1; zmid = zeros(zlen,1);
for i = 1:xlen
  xmid(i,1) = (xedge(i)+xedge(i+1))/2;
end
for i = 1:ylen
  ymid(i,1) = (yedge(i)+yedge(i+1))/2;
end
for i = 1:zlen
  zmid(i,1) = (zedge(i)+zedge(i+1))/2;
end

loc = zeros(size(PCData.lidarDataArray(:,1:3)));
[~,loc(:,1)] = histc(PCData.lidarDataArray(:,1),xedge);
[~,loc(:,2)] = histc(PCData.lidarDataArray(:,2),yedge);
[~,loc(:,3)] = histc(PCData.lidarDataArray(:,3),zedge);

[vCoord,IA,IC] = unique(loc,'rows');
PCData.lidarDataArray(:,9) = IC;
PCData.lidarDataArrayComplete(:,9) = IC;

maxDen = 0;
for i = 1:1:size(vCoord,1)
    selectedIndice = IC==IC(IA(i));
    if(sum(selectedIndice)>maxDen)
        maxDen = sum(selectedIndice);
    end
end

projDataDensified = []; threshPointDesnity = maxDen;
for i = 1:1:size(vCoord,1)
    selectedIndice = IC==IC(IA(i));
    cellData = PCData.lidarDataArray(selectedIndice,:);  
    if(sum(selectedIndice)<3)
        projDataDensified = [projDataDensified; cellData];
    else  
        lackingPointCount = threshPointDesnity - size(cellData,1);
        
        if(lackingPointCount>0)
            randindices = int32(round(rand(lackingPointCount,1)*(size(cellData,1)-1),0))+1;
        else
            randindices = int32(round(rand(threshPointDesnity,1)*(size(cellData,1)-1),0))+1;
        end
        projDataDensified = [projDataDensified; cellData(randindices,:)];
    end
end

%projDataDensified = increasePointDensityToUDValue(projDataDensified);

end

function newLiDARArr = increasePointDensityToUDValue(lidarDataArr)
    numDivisions = 5;
    
    % find layer with max density
    layerPointDensity = zeros(numDivisions,1); 
    mintreeHeight = min(lidarDataArr(:,1));
    maxtreeHeight = max(lidarDataArr(:,1));
    %[mintreeHeight,maxtreeHeight] = deal(num2cell([min(lidarDataArr(:,3)),max(lidarDataArr(:,3))]));    
    divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
    layerPointDensity = 0;
    for i = 1:1:size(divsArr,2)-1
        layerPointDensity(i) = size(lidarDataArr(and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1)),:),1);
    end
    %[maxDensity,maxDensityIndex] = max(layerPointDensity);
    maxDensity = layerPointDensity(4); %round(maxDensity,0);
    constDesnity = 20000;
    newLiDARArr = []; 
    % increase density of remaing layers to the maxDensity
    for i = 1:1:size(divsArr,2)-2
        currentCloud = lidarDataArr(and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1)),:);
        lackingPointCount = constDesnity - size(currentCloud,1) ; 
        randindicesBase = int32(round(rand(constDesnity,1)*(size(constDesnity,1)-1),0))+1;
        if(lackingPointCount>constDesnity)
            randindices = int32(round(rand(lackingPointCount,1)*(size(currentCloud,1)-1),0))+1;
            newLiDARArr = [newLiDARArr; currentCloud(randindicesBase,:); currentCloud(randindices,:)];   
        else
            randindices = int32(round(rand(maxDensity,1)*(size(currentCloud,1)-1),0))+1;
            newLiDARArr = [newLiDARArr; currentCloud(randindicesBase,:); currentCloud(randindices,:)];   
        end
        
    end
    
end
% 
% function newLiDARArr = reducePointHeight(lidarDataArr)
%     numDivisions = 50;
% 
%     mintreeHeight = min(lidarDataArr(:,1)); maxtreeHeight = max(lidarDataArr(:,1));    
%     minL = min(lidarDataArr(:,2)); maxL = max(lidarDataArr(:,2));
%     
%     %maxD = max(lidarDataArr(:,3)); 
%     divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
%     divsArrL = linspace(minL,maxL,numDivisions);
%     
%     newLiDARArr = []; 
%     % increase density of remaing layers to the maxDensity
%     for i = 1:1:size(divsArr,2)-1
%         for j = 1:1:size(divsArrL,2)-1
%             if or( i<(size(divsArr,2)-1), j<(size(divsArrL,2)-1) )                
%                 cond1 = and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1));
%                 cond2 = and(lidarDataArr(:,2) >= divsArrL(j), lidarDataArr(:,2) < divsArrL(j+1));                
%                 currentCloud = lidarDataArr(and(cond1,cond2),:);        
%             else
%                 cond1 = and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) <= divsArr(i+1));
%                 cond2 = and(lidarDataArr(:,2) >= divsArrL(j), lidarDataArr(:,2) <= divsArrL(j+1));                
%                 currentCloud = lidarDataArr(and(cond1,cond2),:);  
%             end
%             oo = currentCloud;
%             oo(:,3) = oo(:,3) - min(oo(:,3)); %+ (maxD - max(oo(:,3)));
%             %oo(:,3) = oo(:,3) + (max(oo(:,3)) - min(oo(:,3)));
%             newLiDARArr = [newLiDARArr; oo ];
%         end
%     end
%     
% end

function newLiDARArr = reducePointHeight(lidarDataArr)
    numDivisions = 50;

    mintreeHeight = min(lidarDataArr(:,1)); maxtreeHeight = max(lidarDataArr(:,1));    
    minL = min(lidarDataArr(:,2)); maxL = max(lidarDataArr(:,2));
    
    medD = median(lidarDataArr(:,3)); 
    divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
    divsArrL = linspace(minL,maxL,numDivisions);
    
    newLiDARArr = []; 
    % increase density of remaing layers to the maxDensity
    for i = 1:1:size(divsArr,2)-1
        for j = 1:1:size(divsArrL,2)-1
            if or( i<(size(divsArr,2)-1), j<(size(divsArrL,2)-1) )                
                cond1 = and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1));
                cond2 = and(lidarDataArr(:,2) >= divsArrL(j), lidarDataArr(:,2) < divsArrL(j+1));                
                currentCloud = lidarDataArr(and(cond1,cond2),:);        
            else
                cond1 = and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) <= divsArr(i+1));
                cond2 = and(lidarDataArr(:,2) >= divsArrL(j), lidarDataArr(:,2) <= divsArrL(j+1));                
                currentCloud = lidarDataArr(and(cond1,cond2),:);  
            end
            oo = currentCloud;            
            diff = (medD - max(oo(:,3)));
            if( and(~isempty(diff), diff>0 ))
                oo(:,3) = oo(:,3) + (medD - max(oo(:,3))); %min(oo(:,3)); % (maxD - max(oo(:,3))); 
            end
                
            %oo(:,3) = oo(:,3) + (max(oo(:,3)) - min(oo(:,3)));
            newLiDARArr = [newLiDARArr; oo ];
        end
    end
    
end

% function newLiDARArr = increasePointDensityToUDValuesss(lidarDataArr,isploton)
%     numDivisions = 50;
%     
%     % find layer with max density
%     layerPointDensity = zeros(numDivisions,1); 
%     mintreeHeight = min(lidarDataArr(:,1));
%     maxtreeHeight = max(lidarDataArr(:,1));
%     
%     minL = min(lidarDataArr(:,2));
%     maxL= max(lidarDataArr(:,2));
%     
%     maxD = max(lidarDataArr(:,3));
%     
%     %[mintreeHeight,maxtreeHeight] = deal(num2cell([min(lidarDataArr(:,3)),max(lidarDataArr(:,3))]));    
%     divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
%     divsArrL = linspace(minL,maxL,numDivisions);
%     
%     newLiDARArr = []; 
%     % increase density of remaing layers to the maxDensity
%     for i = 1:1:size(divsArr,2)-1
%         for j = 1:1:size(divsArrL,2)-1
%             if or( i<(size(divsArr,2)-1), j<(size(divsArrL,2)-1) )                
%                 cond1 = and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1));
%                 cond2 = and(lidarDataArr(:,2) >= divsArrL(j), lidarDataArr(:,2) < divsArrL(j+1));                
%                 currentCloud = lidarDataArr(and(cond1,cond2),:);        
%             else
%                 cond1 = and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) <= divsArr(i+1));
%                 cond2 = and(lidarDataArr(:,2) >= divsArrL(j), lidarDataArr(:,2) <= divsArrL(j+1));                
%                 currentCloud = lidarDataArr(and(cond1,cond2),:);  
%             end
%             oo = currentCloud;
%             oo(:,3) = oo(:,3) + (maxD - max(oo(:,3)));
%             newLiDARArr = [newLiDARArr; oo];
%         end
%     end
%     
% end
% 
% 
% function [projLiDARdata, flippedprojLiDARdata, EV_ratio_thresh] = getProjectedLiDARdata(stData, scaleFactor)
%     
%     plotProjectedData = true; 
%     %plotXScale = 0:1:stData.maxtreeHeight+0.5; 
%     %plotYScale = -10:0.1:10;
%     projLiDARdata = zeros(size(stData.lidarDataArray,1),8);
%     
%     maxXYZ = stData.retMaxXYZ;
%     %retMaxPointIndex = stData.retMaxPointIndexNew;
% 
%     [REFpoint, EV_ratio_thresh]  = getPCAVectors(stData); 
%     %REFpoint(1) = REFpoint(1) - 3;
%     opOffseted = [REFpoint(1:2,:) [10;10] [0;0] [0;0] [stData.lidarDataArray(size(stData.lidarDataArray,1),6)+1;stData.lidarDataArray(size(stData.lidarDataArray,1),6)+1] ];
%     
%     orginalDataXYZ = [stData.lidarDataArray(:,1:3); opOffseted(1,1:3)];
%     refPointXYZ = repmat(opOffseted(1,1:3), size(orginalDataXYZ,1),1);
%     %refPointXYZ = repmat([636049.125552 5129099.038770 57.346000],size(orginalDataXYZ,1),1);
%     
%     dataArrLidAR = [stData.lidarDataArray(:,1:6);opOffseted(1,1:6)];
% 
%     [~, cc]= min(abs(sqrt(sum((orginalDataXYZ.^2-refPointXYZ.^2), 2))));
%     indexExtreme = cc;
%     retVal = [dataArrLidAR(indexExtreme,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexExtreme,3)];
%     
%     indexTreeTop = stData.retMaxXYZIndex;
%     retTreeTop = [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexTreeTop,3)];
%     
%     % Get the plain points for performing calcualtons
%     tempPPa = [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexTreeTop,2) 0];
%     
%     [optimalRefPoint,valleyPointMean] = getOptimalPlaneRefPoint(stData);
%     if(~isempty(optimalRefPoint))
%         %EV_ratio_thresh = median(valleyPointMean);
%         if(optimalRefPoint(1)<0)
%             optimalRefPoint = -optimalRefPoint;
%         end
%     end
%     
%     if(isempty(optimalRefPoint))
%         tempPPb = [dataArrLidAR(indexExtreme,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexExtreme,3)];  
%     else
%         tempPPb = [optimalRefPoint(1,1) optimalRefPoint(1,2) dataArrLidAR(indexExtreme,3)];
%     end
%     %tempPPb = [dataArrLidAR(indexExtreme,1) dataArrLidAR(indexExtreme,2) dataArrLidAR(indexExtreme,3)];
%     
%     tempPPc = [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexTreeTop,2) dataArrLidAR(indexTreeTop,3)];
%     
%     %ppa = [retTreeTop(1) retTreeTop(2) 0];
%     %ppb = [stData.lidarDataArray(indexExtreme,1) stData.lidarDataArray(indexExtreme,2) stData.lidarDataArray(indexExtreme,3)];
%     %ppc = retTreeTop;
%     sz = size(dataArrLidAR,1);
%     for i = 1:1:sz          
%         inputPoint = dataArrLidAR(i,1:3);            
%         inputPointIndex = dataArrLidAR(i,6);  
%         externPoint = inputPoint;        
%         position = getPointDirection(tempPPa, tempPPb, tempPPc, externPoint);
%         [sectorLength, disttocentre] = getSectorLength(tempPPb, inputPoint, [dataArrLidAR(indexTreeTop,1) dataArrLidAR(indexTreeTop,2)], [maxXYZ(1) maxXYZ(2)], scaleFactor); %% disttoCentre = Y, sectorlength = X,  inputPointZ = Z  
%         projLiDARdata(i,1:8) = [inputPoint(3) position*sectorLength  disttocentre 0 inputPoint inputPointIndex];    
%     end
%     
%     projLiDARdata = real(projLiDARdata);
%     
%     if(and(NC.ISPLOTON,plotProjectedData))
%         f2 = figure;  
%         set(f2, 'Position', [10 60 600 420])           
%         subplot(1,2,1);
%         hold on; 
%         uu = scatter3(projLiDARdata(:,1), projLiDARdata(:,2), projLiDARdata(:,3),300, projLiDARdata(:,3),'.');
%         xlabel('Tree Height'); ylabel('Distance to the refernce cyclinder surface'); zlabel('Distance to Stem');
%         camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-90,90); box on;
%         axis([min(projLiDARdata(:,1)) max(projLiDARdata(:,1)) min(projLiDARdata(:,2)) max(projLiDARdata(:,2)) min(projLiDARdata(:,3)) max(projLiDARdata(:,3))]);
%         title('Projected LiDAR Data');
%         T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
%         x = [0 85 110 255];
%         map = interp1(x/255,T,linspace(0,1,100));
%         colormap(map)
%         C = uu.CData;       
%     end
% 
%     % Get the points for making the 3D plane
%     ppa = [0 0 projLiDARdata(indexExtreme,3)];
%     ppb = [projLiDARdata(indexTreeTop,1) 0 projLiDARdata(indexTreeTop,3)];
%     ppc = [projLiDARdata(indexExtreme,1) 0 projLiDARdata(indexExtreme,3)];    
%     plainPointsAll = [ppa;ppb;ppc;ppa];
%     
%     if(and(NC.ISPLOTON,plotProjectedData))
%         hold on;
%         subplot(1,2,1);
%         plot3(plainPointsAll(:,1), plainPointsAll(:,2) ,plainPointsAll(:,3), '.', 'MarkerSize', 40,'Color', [1 0 0]);
%         hold on;    
%         plotPlaneFrom3Points(ppa, ppb, ppc, false);
%         text(ppa(1),ppa(2),ppa(3),'Projection of tree top on XY plane','HorizontalAlignment','right');
%         text(ppc(1),ppc(2),ppc(3),'reference point','HorizontalAlignment','right');
%         text(ppb(1),ppb(2),ppb(3),'Tree top', 'HorizontalAlignment','right');
%     end
%     projLiDARdata(:,4) = [1:1:size(stData.lidarDataArray,1)+1]';
%     
%     % remove reference point added
%     projLiDARdata = projLiDARdata(1:end-1,:);
%     
%     % flattern cloud
% % %    mintreeHeight = min(projLiDARdata(:,1));  maxtreeHeight = max(projLiDARdata(:,1));
% % %    divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);    
% %     maxD = max(projLiDARdata(:,3));    
% % %    newLiDARArr = []; 
% %     for i = 1:1:size(divsArr,2)-1
% %         cond1 = and(projLiDARdata(:,1) >= divsArr(i), projLiDARdata(:,1) < divsArr(i+1));          
% %         currentCloud = projLiDARdata(and(cond1,cond2),:);
% %         currentCloud(:,3) = currentCloud(:,3)+(maxD-max(currentCloud(:,3)));
% %         newLiDARArr = [newLiDARArr; currentCloud];
% %     end
% %     projLiDARdata = newLiDARArr;
%     %projLiDARdata = reducePointHeight(projLiDARdata);
% 
%     
%     % new analysis
%     flippedprojLiDARdata = projLiDARdata;
%     %flippedprojLiDARdata = increasePointDensityToUDValuesss(projLiDARdata);
%     %flippedprojLiDARdata(:,3) = max(flippedprojLiDARdata(:,3)) - flippedprojLiDARdata(:,3);  
%     
%     if(and(NC.ISPLOTON,plotProjectedData))    
%         subplot(1,2,2);hold on;
%         uu = scatter3(flippedprojLiDARdata(:,1), flippedprojLiDARdata(:,2), flippedprojLiDARdata(:,3),300, flippedprojLiDARdata(:,3),'.');
%         xlabel('Tree Height'); ylabel('Distance to the refernce cyclinder surface'); zlabel('Distance to Stem');
%         camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-90,90); box on;
%         axis([min(projLiDARdata(:,1)) max(projLiDARdata(:,1)) min(projLiDARdata(:,2)) max(projLiDARdata(:,2)) min(projLiDARdata(:,3)) max(projLiDARdata(:,3))]);
%         title('Projected LiDAR Data');
%         T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
%         x = [0 85 110 255];
%         map = interp1(x/255,T,linspace(0,1,100));
%         colormap(map);
%         C = uu.CData;
%     end
%     
%     if(and(NC.ISPLOTON,flippedprojLiDARdata))
%         subplot(1,2,2);
%         plot3(plainPointsAll(:,1), plainPointsAll(:,2) ,plainPointsAll(:,3), '.', 'MarkerSize', 40,'Color', [1 0 0]);
%         hold on;    
%         plotPlaneFrom3Points(ppa, ppb, ppc, false);
%         text(ppa(1),ppa(2),ppa(3),'Projection of tree top on XY plane','HorizontalAlignment','right');
%         text(ppc(1),ppc(2),ppc(3),'reference point','HorizontalAlignment','right');
%         text(ppb(1),ppb(2),ppb(3),'Tree top', 'HorizontalAlignment','right');
%     end
%     
% end
% 
% function projDataDensified = increasePointDensityToUDValuesss(PJData)
%    PJData = projDataPreProcessing(PJData); 
%    VDim = getVoxelDimension(PJData, 0.5);
%    projDataDensified = getVoxelDensityAndNebrCellIndicesNewsss(PJData, VDim.rX, VDim.rY, VDim.rZ);
% end
% 
% function projDataDensified = getVoxelDensityAndNebrCellIndicesNewsss(PCData, xedge, yedge, zedge)
% 
% xlen = length(xedge)-1; xmid = zeros(xlen,1);
% ylen = length(yedge)-1; ymid = zeros(ylen,1);
% zlen = length(zedge)-1; zmid = zeros(zlen,1);
% for i = 1:xlen
%   xmid(i,1) = (xedge(i)+xedge(i+1))/2;
% end
% for i = 1:ylen
%   ymid(i,1) = (yedge(i)+yedge(i+1))/2;
% end
% for i = 1:zlen
%   zmid(i,1) = (zedge(i)+zedge(i+1))/2;
% end
% 
% loc = zeros(size(PCData.lidarDataArray(:,1:3)));
% [~,loc(:,1)] = histc(PCData.lidarDataArray(:,1),xedge);
% [~,loc(:,2)] = histc(PCData.lidarDataArray(:,2),yedge);
% [~,loc(:,3)] = histc(PCData.lidarDataArray(:,3),zedge);
% 
% [vCoord,IA,IC] = unique(loc,'rows');
% PCData.lidarDataArray(:,9) = IC;
% PCData.lidarDataArrayComplete(:,9) = IC;
% 
% maxDen = 0;
% for i = 1:1:size(vCoord,1)
%     selectedIndice = IC==IC(IA(i));
%     if(sum(selectedIndice)>maxDen)
%         maxDen = sum(selectedIndice);
%     end
% end
% 
% projDataDensified = []; threshPointDesnity = maxDen;
% for i = 1:1:size(vCoord,1)
%     selectedIndice = IC==IC(IA(i));
%     cellData = PCData.lidarDataArray(selectedIndice,:);  
%     if(sum(selectedIndice)<3)
%         projDataDensified = [projDataDensified; cellData];
%     else  
%         lackingPointCount = threshPointDesnity - size(cellData,1);
%         
%         if(lackingPointCount>0)
%             randindices = int32(round(rand(lackingPointCount,1)*(size(cellData,1)-1),0))+1;
%         else
%             randindices = int32(round(rand(threshPointDesnity,1)*(size(cellData,1)-1),0))+1;
%         end
%         projDataDensified = [projDataDensified; cellData(randindices,:)];
%     end
% end
% 
% %projDataDensified = increasePointDensityToUDValue(projDataDensified);
% 
% end
% 
% function newLiDARArr = increasePointDensityToUDValue(lidarDataArr)
%     numDivisions = 5;
%     
%     % find layer with max density
%     layerPointDensity = zeros(numDivisions,1); 
%     mintreeHeight = min(lidarDataArr(:,1));
%     maxtreeHeight = max(lidarDataArr(:,1));
%     %[mintreeHeight,maxtreeHeight] = deal(num2cell([min(lidarDataArr(:,3)),max(lidarDataArr(:,3))]));    
%     divsArr = linspace(mintreeHeight,maxtreeHeight,numDivisions);
%     layerPointDensity = 0;
%     for i = 1:1:size(divsArr,2)-1
%         layerPointDensity(i) = size(lidarDataArr(and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1)),:),1);
%     end
%     %[maxDensity,maxDensityIndex] = max(layerPointDensity);
%     maxDensity = layerPointDensity(4); %round(maxDensity,0);
%     constDesnity = 20000;
%     newLiDARArr = []; 
%     % increase density of remaing layers to the maxDensity
%     for i = 1:1:size(divsArr,2)-2
%         currentCloud = lidarDataArr(and(lidarDataArr(:,1) >= divsArr(i), lidarDataArr(:,1) < divsArr(i+1)),:);
%         lackingPointCount = constDesnity - size(currentCloud,1) ; 
%         randindicesBase = int32(round(rand(constDesnity,1)*(size(constDesnity,1)-1),0))+1;
%         if(lackingPointCount>constDesnity)
%             randindices = int32(round(rand(lackingPointCount,1)*(size(currentCloud,1)-1),0))+1;
%             newLiDARArr = [newLiDARArr; currentCloud(randindicesBase,:); currentCloud(randindices,:)];   
%         else
%             randindices = int32(round(rand(maxDensity,1)*(size(currentCloud,1)-1),0))+1;
%             newLiDARArr = [newLiDARArr; currentCloud(randindicesBase,:); currentCloud(randindices,:)];   
%         end
%         
%     end
%     
% end
% 




% function newCloud = raisePointsinCSM(currentCloud)
%     newCloud = [];
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
%     medDHeight = median(currentCloud(:,3));
%     
%     cnt1 = 1;
%     for i = MMS.minX+halfXGridSize:gXStep:MMS.maxX-halfXGridSize
%         cnt2 = 1;
%         for j = MMS.minY+halfYGridSize:gYStep:MMS.maxY-halfYGridSize
%             cond1 = and(currentCloud(:,1)>i-halfXGridSize, currentCloud(:,1)<=i+halfXGridSize);
%             cond2 = and(currentCloud(:,2)>j-halfYGridSize, currentCloud(:,2)<=j+halfYGridSize);
%             indGridPoints = and(cond1,cond2);            
%             
%             localMaxheight = max(currentCloud(indGridPoints,3));            
%             diffHeight = (medDHeight - localMaxheight);
%             
%             if(~isempty(diffHeight))
%                 if(diffHeight>0)           
%                     newCloud = [newCloud; [currentCloud(indGridPoints,1:2) currentCloud(indGridPoints,3)+diffHeight ] ];
%                 else
%                     newCloud = [newCloud; currentCloud(indGridPoints,1:3)];
%                 end
%             end
%             cnt2 = cnt2+1;
%         end
%         cnt1 = cnt1+1;
%     end   
% end
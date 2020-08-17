function [csmcdmNormalized,dsmImageSegmentedArr, colorArr] =  getCSMCDM(threshold, slabCord, VOXDiv, filterOrder, plotCSM)
    inputParams.DimGrid1 = size(VOXDiv.rX,2); inputParams.DimGrid2 = size(VOXDiv.rY,2);
    inputParams.gDim1Step = VOXDiv.xStep; inputParams.gDim2Step = VOXDiv.yStep;
    inputParams.halfDim1GridSize = VOXDiv.xStep*1; inputParams.halfDim2GridSize = VOXDiv.yStep*5;
    inputParams.filterOrder = filterOrder;
    inputParams.threshold = threshold;
    inputParams.plotCSM = plotCSM;
    inputParams.slabCord = slabCord;
    MMS = getMaxMins(slabCord);
    inputParams.minDim1 = MMS.XMin; inputParams.maxDim1 = MMS.XMax;
    inputParams.minDim2 = MMS.YMin; inputParams.maxDim2 = MMS.YMax;
    inputParams.plotSurfCSM = false;
    [csmcdmNormalized,dsmImageSegmentedArr,colorArr] = getCSMCDMfn(inputParams);
end

function [csmcdm,dsmImageSegmentedArr,colorArr] = getCSMCDMfn(inP)
    % Initializing csm and cdm arrays
    csm = zeros(inP.DimGrid1,inP.DimGrid2);
    cdm = zeros(inP.DimGrid1,inP.DimGrid2);
    % Generating csm and cdm
    cnt1 = 1;
    for i = inP.minDim1+inP.halfDim1GridSize:inP.gDim1Step:inP.maxDim1-inP.halfDim1GridSize
        cnt2 = 1;
        for j = inP.minDim2+inP.halfDim2GridSize:inP.gDim2Step:inP.maxDim2-inP.halfDim2GridSize
            cond1 = and(inP.slabCord(:,1)>i-(inP.halfDim1GridSize), inP.slabCord(:,1)<=i+(inP.halfDim1GridSize));
            cond2 = and(inP.slabCord(:,2)>j-(inP.halfDim2GridSize), inP.slabCord(:,2)<=j+(inP.halfDim2GridSize));
            indGridPoints = find(and(cond1,cond2));
            if(size(indGridPoints,1)>0)
                csm(cnt1,cnt2) = mean(inP.slabCord(indGridPoints,3));
                cdm(cnt1,cnt2) = length(indGridPoints)*(max(inP.slabCord(indGridPoints,3)) - min(inP.slabCord(indGridPoints,3)));
                %imagesc(cdm); pause(0.02); hold on;
            end
            cnt2 = cnt2+1;
        end
        cnt1 = cnt1+1;
    end

   % Intepolate CSM and CSM where values are missing
   csm = interploteImage(csm, 'linear');
   cdm = interploteImage(cdm, 'linear'); % Mean Median linear
  
   %filterAnalysis(csm) % for test
   inP.filterOrder = 2;
   csm = imgaussfilt(csm,inP.filterOrder);
   cdm = imgaussfilt(cdm,inP.filterOrder);
   
   % remove noisy points due to isolated branches
   se = strel('disk',5);
   cdm = imerode(cdm,se);
   cdm = imdilate(cdm,se);
   cdm = imclose(cdm,se);
   
   % smoothen CSM using vertical-rectangle filter (as trees grows vertically)
   csm = imgaussfilt(csm,[round(size(csm,1)*0.1,0) 3]);
   cdm = imgaussfilt(cdm,[round(size(cdm,1)*0.1,0) 3]);
   cdm=mat2gray(cdm); csm=mat2gray(csm);
   
   % merge csm and cdm (multiplicatin instead of averging improved sub-tree visibility)
   csmcdm = flipud(mat2gray(cdm.*csm)); %max(csm,cdm) %  cdm .*
   %[counts,~] = imhist(csmcdm,256); thresholdfromcsmcdm = otsuthresh(counts);

   % repeat the kmeans clustering 3 times to avoid local minima
   KMInput = [csm(:) cdm(:)];
   KMInput(isnan(KMInput)) = 0;
   noOfClasses = 3;
   [cluster_idx, cluster_center] = kmeans(KMInput,noOfClasses,'distance','sqEuclidean','Replicates',3);
   cluster_idx = reshape(cluster_idx,size(csm,1),size(csm,2));   
   cvalsArr = unique(cluster_idx);
   SumArr = [];
   for cid = 1:1:size(cvalsArr)
      indxofclass = cluster_idx(:)==cvalsArr(cid);
      tempcsmcdm = flipud(csmcdm);
      SumArr = [SumArr; [(sum(tempcsmcdm(indxofclass))/sum(indxofclass)) cvalsArr(cid)]];
   end
   [maxSum, maxSumIndex] = max(SumArr(:,1));
   % take class with max values in CSMCDM pixels, and merge remaining classes
   ClassofInterest = SumArr(SumArr(:,2)==maxSumIndex,2);
   ClassofDisinterest = setdiff(cvalsArr,ClassofInterest);
   for j = 1:1:length(ClassofDisinterest)
        cluster_idx(cluster_idx(:)==ClassofDisinterest(j)) = 0;
   end
   % plot segmented map
   segmentedcsmcdm = flipud(cluster_idx);
      
   % Get subdominant tree boundary from CSM
   plotOn = true;
   [dsmImageSegmentedArr, colorArr] = getSubTreeBoundingBox(flipud(segmentedcsmcdm), csmcdm, inP.threshold, plotOn);
   dsmImageSegmented = zeros(size(dsmImageSegmentedArr{1}));
   for idsmCnt=1:1:size(dsmImageSegmentedArr,2)
        dsmImageSegmented = dsmImageSegmented + dsmImageSegmentedArr{idsmCnt};
   end
   
   % plot the tree boundary
   if(and(NC.ISPLOTON, inP.plotCSM))
        f6=figure('name','CSM Plots');
        set(f6, 'Position', [690 60 600 420]);
        imagesc(csmcdm); hold on;
        %dsmImageSegmented = flip(dsmImageSegmented);
        imagesc(mat2gray(dsmImageSegmented)); alpha(0.5);
        colormap(jet);
        [heightY,lengthX] = size(csm);
        dd = lengthX/2;
        xTickArr = -dd-(5-mod(dd,5)):5:dd+(5-mod(dd,5));
        yTickArr = 0:5:heightY+(5-mod(heightY,5));
        set(gca,'XTick', xTickArr+dd+(5-mod(dd,5)));
        set(gca,'XTickLabel', num2cell(ceil(xTickArr*inP.gDim1Step)));
        set(gca,'YTick', yTickArr);
        set(gca,'YTickLabel', num2cell( floor( flip(mat2gray(yTickArr)*inP.maxDim2)) ));
        set(findall(gcf,'type','axes'),'fontsize',32);
        set(findall(gcf,'type','text'),'fontSize',32);
        ylabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 36);
        xlabel('Distance to Reference Point','Fontname', 'Times New Roman' ,'FontSize', 36);
        title('Crown Surface Model');
        %axis equal;
   end

end

% function [csmNormalized,dsmImageSegmentedArr, colorArr] =  getCSMCDM(threshold, slabCord, FProjDat, VOXDiv, filterOrder, plotCSM)
% 
%     colorArr =[];
%     MMS = getMaxMins(slabCord);
%     VOXEL_XYZ_DIVS = getVoxelDimension(MMS, 0.25);
% 
%     gXStep = abs(MMS.maxX - MMS.minX)/(VOXEL_XYZ_DIVS.xDiv);
%     gYStep = abs(MMS.maxY - MMS.minY)/(VOXEL_XYZ_DIVS.yDiv);   
%     gZStep = abs(MMS.maxZ - MMS.minZ)/(VOXEL_XYZ_DIVS.zDiv);   
% 
%     halfXGridSize  = gXStep/2;       
%     halfYGridSize  = gYStep/2;
%     halfZGridSize  = gZStep/2; 
% 
%     xDimGrid = size(MMS.minX+halfXGridSize : gXStep : MMS.maxX-halfXGridSize,2);
%     yDimGrid = size(MMS.minY+halfYGridSize : gYStep : MMS.maxY-halfYGridSize,2);        
%     zDimGrid = size(MMS.minZ+halfZGridSize : gZStep : MMS.maxZ-halfZGridSize,2);    
% 
%     inputParams.DimGrid1 = xDimGrid;    inputParams.DimGrid2 = yDimGrid;
%     inputParams.gDim1Step = gXStep;    inputParams.gDim2Step = gYStep;
%     inputParams.minDim1 = MMS.XMin;    inputParams.maxDim1 = MMS.XMax;
%     inputParams.minDim2 = MMS.YMin;    inputParams.maxDim2 = MMS.YMax;
%     inputParams.halfDim1GridSize = halfXGridSize;    inputParams.halfDim2GridSize = halfYGridSize;
%     inputParams.slabCord = slabCord;    inputParams.filterOrder = filterOrder;
%     inputParams.threshold = threshold;    inputParams.plotCSM = plotCSM;
%     inputParams.FProjDat = FProjDat;
%     
%     [csmNormalized,dsmImageSegmentedArr] = getCSMCDMfn(inputParams);
% end

% function [csmNormalized,dsmImageSegmentedArr] = getCSMCDMfn(inputParams)
% figure;
%     DimGrid1 = inputParams.DimGrid1;
%     DimGrid2 = inputParams.DimGrid2;
%     gDim1Step = inputParams.gDim1Step;
%     gDim2Step = inputParams.gDim2Step;
%     minDim1 = inputParams.minDim1;
%     maxDim1 = inputParams.maxDim1;
%     minDim2 = inputParams.minDim2;
%     maxDim2 = inputParams.maxDim2;
%     halfDim1GridSize = inputParams.halfDim1GridSize;
%     halfDim2GridSize = inputParams.halfDim2GridSize;
%     slabCord = inputParams.slabCord;
%     filterOrder = inputParams.filterOrder;
%     threshold = inputParams.threshold;
%     plotCSM = inputParams.plotCSM;
%     FProjDat = inputParams.FProjDat;
% 
%     % Initializing csm and cdm arrays
%     csm = zeros(DimGrid1,DimGrid2);
%     %csmNew = zeros(DimGrid1,DimGrid2);
%     cdm = zeros(DimGrid1,DimGrid2);
% 
%     % Generating csm and cdm
%     cnt1 = 1;
%     for i = minDim1+halfDim1GridSize:gDim1Step:maxDim1-halfDim1GridSize
%         cnt2 = 1;
%         for j = minDim2+halfDim2GridSize:gDim2Step:maxDim2-halfDim2GridSize               
%             cond1 = and(slabCord(:,1)>i-halfDim1GridSize, slabCord(:,1)<=i+halfDim1GridSize);
%             cond2 = and(slabCord(:,2)>j-halfDim2GridSize, slabCord(:,2)<=j+halfDim2GridSize);
%             cond3 = true;%slabCord(:,3)>(max(slabCord(:,3))*0.5);
%             indGridPoints = find(and(and(cond1,cond2),cond3));
% 
%             if(size(indGridPoints,1)>0)  
%                 csm(cnt1,cnt2) = mean(slabCord(indGridPoints,3));
%                 %csmNew(cnt1,cnt2) = mean(FProjDat(indGridPoints,3));
%                 %aa = rangesearch(slabCord(:,1:2),slabCord(indGridPoints,1:2),4);                
%                 %cdm(cnt1,cnt2) = size(aa{1},2);
%                 cdm(cnt1,cnt2) = length(indGridPoints)*(max(slabCord(indGridPoints,3)) - min(slabCord(indGridPoints,3)));
%             end
%             cnt2 = cnt2+1;
%         end
%         cnt1 = cnt1+1;
%     end
%     
%    csm = ((csm-min(csm(:)))/( max(csm(:)) - min(csm(:)) ))*255;
%    %csmNew = ((csmNew-min(csmNew(:)))/( max(csmNew(:)) - min(csmNew(:)) ))*255;
%    cdm = ((cdm-min(cdm(:)))/( max(cdm(:)) - min(cdm(:)) ))*255;
% 
%    % Intepolate CSM where values are missing
%    idx = find(csm==0);        
%    [r,c] = ind2sub(size(csm),idx);      
%    for j = 1:1:size(r,1)
%       avgVal = getRC(r(j),c(j),csm);
%       csm(r(j),c(j)) = avgVal;
%    end
% %    % Intepolate csmNew where values are missing
% %    idx = find(csmNew==0);        
% %    [r,c] = ind2sub(size(csmNew),idx);      
% %    for j = 1:1:size(r,1)
% %       avgVal = getRC(r(j),c(j),csmNew);
% %       csmNew(r(j),c(j)) = avgVal;
% %    end
%    % Intepolate CDM where values are missing
%    idx = find(cdm==0);        
%    [r,c] = ind2sub(size(cdm),idx);      
%    for j = 1:1:size(r,1)
%       avgVal = getRC(r(j),c(j),cdm);
%       cdm(r(j),c(j)) = avgVal;
%    end
%    
%    % Gaussian smotthening CSM and CDM
%    csm = imgaussfilt(fliplr(flip(csm)),filterOrder);
%    %csmNew = imgaussfilt(fliplr(flip(csmNew)),filterOrder);
%    cdm = imgaussfilt(fliplr(flip(cdm)),filterOrder);
%    % merge information in csm and csm (pixel wise mean)
%    %csm =cdm;
%    %csmNew = (csmNew); %(csm + cdm)/2;
%    
%    %csmNormalized11=mat2gray(csmNew);
%    
%    csm = imgaussfilt(csm,[round(size(csm,1)*0.05,0) 1]);
%    cdm = imgaussfilt(cdm,[round(size(cdm,1)*0.1,0) 3]);
%    
%    subplot(2,3,1); surf(cdm); title('cdm');
%    subplot(2,3,2); surf(csm); title('csm');
%    %subplot(2,3,3); surf(csmNew); title('csmNew');
%    cdm=mat2gray(cdm);
%    csm=mat2gray(csm);
%    %csmNew=mat2gray(csmNew);
%    subplot(2,3,5); surf(min(csm,cdm)); title('csmCDM');
%    %subplot(2,3,6); surf(min(csmNew,cdm)); title('csmNewCDM');
%    
%    csmNormalized=mat2gray(min(csm,cdm));
%    % get approximate crown extend
%    plotOn = false;
%   % [counts,x] = imhist(csmNormalized,16); 
%    % T = otsuthresh(counts);
%    
%    [dsmImageSegmentedArr, colorArr] = getSubTreeBoundingBox(csmNormalized, threshold, plotOn);
%    dsmImageSegmented = zeros(size(csm));
%    for idsmCnt=1:1:size(dsmImageSegmentedArr,2)
%         dsmImageSegmented = dsmImageSegmented + dsmImageSegmentedArr{idsmCnt};
%    end
%    
% %    [dsmImageSegmentedArr11, colorArr] = getSubTreeBoundingBox(csmNormalized11, threshold, plotOn);
% %    for idsmCnt=1:1:size(dsmImageSegmentedArr11,2)
% %         dsmImageSegmented = dsmImageSegmented + dsmImageSegmentedArr11{idsmCnt};
% %    end
% 
%    if(and(NC.ISPLOTON, plotCSM))
%         f6=figure('name','CSM Plots');                       
%         set(f6, 'Position', [690 60 600 420])
%         imagesc(csmNormalized); hold on;
%         %dsmImageSegmented = flip(dsmImageSegmented);
%         imagesc(mat2gray(dsmImageSegmented)); alpha(0.5);        
%         colormap(jet);        
%         [heightY,lengthX] = size(csm);
%         dd = lengthX/2;
%         xTickArr = -dd-(5-mod(dd,5)):5:dd+(5-mod(dd,5));
%         yTickArr = 0:5:heightY+(5-mod(heightY,5));
%         set(gca,'XTick', xTickArr+dd+(5-mod(dd,5)));
%         set(gca,'XTickLabel', num2cell(ceil(xTickArr*gDim1Step)));
%         set(gca,'YTick', yTickArr);
%         set(gca,'YTickLabel', num2cell( floor( flip(mat2gray(yTickArr)*maxDim2)) ));        
%         set(findall(gcf,'type','axes'),'fontsize',32)
%         set(findall(gcf,'type','text'),'fontSize',32)
%         ylabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 36); 
%         xlabel('Distance to Reference Point','Fontname', 'Times New Roman' ,'FontSize', 36);
%         title('Crown Surface Model');
%         %axis equal;
%    end
% 
% end
   %csm = ((csm-min(csm(:)))/( max(csm(:)) - min(csm(:)) ))*255;
   %cdm = ((cdm-min(cdm(:)))/( max(cdm(:)) - min(cdm(:)) ))*255;

   
      
%    % Surface plot of CSMCDM
%    if(and(NC.ISPLOTON, inP.plotSurfCSM))
%        figure;
%        subplot(1,3,1); surf(cdm); title('cdm');
%        subplot(1,3,2); surf(csm); title('csm');
%        subplot(1,3,3); surf(csmcdm); title('csmCDM');
%    end
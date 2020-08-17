function [optimalRefPoint, valleyPointMean] = getOptimalPlaneRefPoint(stData)    
    %subplot(1,2,2);   
    
    %plot3(stData.lidarDataArray(:,1),stData.lidarDataArray(:,2),stData.lidarDataArray(:,3),'*');
    %hold on; axis equal;
    k = boundary(stData.lidarDataArray(:,1),stData.lidarDataArray(:,2),0.5);
    boundPoints = [stData.lidarDataArray(k,1) stData.lidarDataArray(k,2)];
    %plot(boundPoints(:,1),boundPoints(:,2),'-*'); hold on;
    
    idx = stData.retMaxXYZIndex;    
    maxXYZ = stData.lidarDataArray(idx,1:2); 
    
    %plot(maxXYZ(1),maxXYZ(2),'.','MarkerSize',35); hold on; axis square;
    distValArr = pdist2(boundPoints, maxXYZ);
    
    distValArrMAvged = movmean(distValArr,20);
    XArr = [1:1:size(distValArrMAvged,1)];
    %plot(XArr,distValArr,'-*'); hold on; axis square;
    
    valleyIndex = peakfinder(distValArrMAvged,0,100,-1, false);

    %plot(XArr(peakIndex),valleyPoints,'*'); hold on; axis square;
    
    optimalRefPoint = sortrows([boundPoints(valleyIndex,:) distValArr(valleyIndex,:)],3);
    
    if(~isempty(optimalRefPoint))
        valleyPointMean = optimalRefPoint(1,3)/stData.treeWidth;
    else
        valleyPointMean = 0.5;
    end
    
    %plot(optimalRefPoint(1,1), optimalRefPoint(1,2), '*'); hold on; axis square;
    %axis([0 120 0 15]);
end
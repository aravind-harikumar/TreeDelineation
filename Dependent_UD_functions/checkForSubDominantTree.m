function retVal = checkForSubDominantTree(pdData,treeCnt)
    projectedData = pdData.lidarDataArray;
    VOXEL_XYZ_DIVS = getVoxelDimension(pdData, 1);
    %histogram(projectedData(:,1),10);
    retVal =0;
    [hst,opHist] = histc( projectedData(:,1),VOXEL_XYZ_DIVS.rX);
    histVal = [];    
    oop = unique(opHist);    
    for i = 1:1:size(oop,1);
        vv = opHist(:)==oop(i);
        histVal = [histVal; max(projectedData(vv,3))];
    end

    
    histVal = histVal - min(histVal);
    
    w = gausswin(5);
    histVal = filter(w,1,histVal);
    %histVal = movmean(histVal,[0 5]);
    histVal = (histVal-min(histVal(:)))/(max(histVal(:))-min(histVal(:)));
    %histVal = diff(histVal); histVal = histVal - min(histVal);
    %histVal = (histVal-min(histVal(:)))/(max(histVal(:))-min(histVal(:)));
    ValleyIndex = peakfinder(histVal,0.05,0,-1, true);
    ValleyXYOrig = [ValleyIndex histVal(ValleyIndex)];

    if(true)
        subplot(10,10,treeCnt);
        plot([1:1:size(histVal,1)]',histVal,'-'); hold on;
        if(~isempty(ValleyXYOrig))
            plot(ValleyXYOrig(:,1),ValleyXYOrig(:,2),'.b', 'MarkerSize', 10);
            stem(ValleyXYOrig(:,1),ValleyXYOrig(:,2));
        end
        axis([0 size(oop,1) 0 1]);
    end
end

%xx = 1:0.5:size(histVal,1);
%histVal = spline([1:1:size(histVal,1)]',histVal,xx);
%histVal = movmean(histVal,5); % x,5,'Endpoints','fill'
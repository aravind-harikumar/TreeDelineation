function retPointIndx = getProminalPoints(startPoint, otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)
    [dist,Indx] = pdist2(otherPoints,startPoint,'euclidean','Smallest',NumNearestPoints);
    retPointIndx = [retPointIndx;Indx];
    maxIter = maxIter-1;
    for i = 2:1:2
        if(dist(i) < distTheshold) && maxIter > 0
            retPointIndx = [retPointIndx;getProminalPoints(otherPoints(retPointIndx(size(retPointIndx,1)),:), otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)];
        else
            retPointIndx = [];
        end
    end  
end
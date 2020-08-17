
function retDistWeights = calculateDistWeights(lidarDataArray, neighbourhoodSize)
   [nearPoint,~] = pdist2(lidarDataArray(:,(1:3)),lidarDataArray(:,(1:3)),'euclidean','Smallest',100);  
    nearPoint(nearPoint>=neighbourhoodSize) = 0; % make dist 0 for points > threshold  neighbourhoodDist  
    nearPoint = nearPoint(2:end,:); % remove first line because it is distance to point itself.
    neigDensityArr = sum(nearPoint~=0,1)'; % find non zero elenents in the array for every column
    retDistWeights = neigDensityArr/norm(neigDensityArr); % normalization
end
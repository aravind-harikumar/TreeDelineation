function [retMaxXYZ, maxpointIndex] = findMaxHeightXYNew(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
     maxpointIndex = lidarDataArray(indx,6);
end

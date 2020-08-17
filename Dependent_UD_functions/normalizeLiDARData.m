function lidarDataArray = normalizeLiDARData(lidarDataArray)
    midxPoint = (max(lidarDataArray(:,1))- min(lidarDataArray(:,1)))/2;
    midyPoint = (max(lidarDataArray(:,2))- min(lidarDataArray(:,2)))/2;
    
    lidarDataArray(:,1) = lidarDataArray(:,1)- min(lidarDataArray(:,1));
    lidarDataArray(:,2) = lidarDataArray(:,2)- min(lidarDataArray(:,2));
    lidarDataArray(:,3) = lidarDataArray(:,3)- min(lidarDataArray(:,3));

    lidarDataArray(:,1) = lidarDataArray(:,1)- midxPoint;
    lidarDataArray(:,2) = lidarDataArray(:,2)- midyPoint;
end
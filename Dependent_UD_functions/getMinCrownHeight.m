function crownHt = getMinCrownHeight(lidarDataArray)
   maxTreeHeight = max(lidarDataArray(:,3));
   crownHt = min(lidarDataArray(:,3));
   for i = maxTreeHeight:-1:1
       templidararray = lidarDataArray(find(and(lidarDataArray(:,3)<i,(lidarDataArray(:,3)>i-1))),1:3);
        if(size(templidararray,1)<10)
            crownHt = i;
            break;
        end
   end
end
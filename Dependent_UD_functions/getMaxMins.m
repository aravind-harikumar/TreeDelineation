function pr = getMaxMins(slabCord)
    pr.maxX = max(slabCord(:,1));
    pr.minX = min(slabCord(:,1));
    pr.maxY = max(slabCord(:,2));
    pr.minY = min(slabCord(:,2));
    pr.maxZ = max(slabCord(:,3));
    pr.minZ = min(slabCord(:,3));
    
    pr.XMax = max(slabCord(:,1));
    pr.XMin = min(slabCord(:,1));
    pr.YMax = max(slabCord(:,2));
    pr.YMin = min(slabCord(:,2));
    pr.ZMax = max(slabCord(:,3));
    pr.ZMin = min(slabCord(:,3));
end

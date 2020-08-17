function idx = incircle(lidarDataArray, point, radius)
    [idx,~] = rangesearch(lidarDataArray,point,radius);
    idx = idx{1};
end
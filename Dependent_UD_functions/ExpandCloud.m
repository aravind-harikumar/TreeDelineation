function branchPointCluster = ExpandCloud(branchPointCluster, scale)
    tmp = branchPointCluster;
    branchPointCluster = branchPointCluster*scale;
    xdiff = (max(branchPointCluster(:,1)) - max(tmp(:,1))) - (max(branchPointCluster(:,1))-min(branchPointCluster(:,1)))*0.5 + (max(tmp(:,1))-min(tmp(:,1)))*0.5;
    ydiff = (max(branchPointCluster(:,2)) - max(tmp(:,2))) - (max(branchPointCluster(:,2))-min(branchPointCluster(:,2)))*0.5 + (max(tmp(:,2))-min(tmp(:,2)))*0.5;
    zdiff = (max(branchPointCluster(:,3)) - max(tmp(:,3))) - (max(branchPointCluster(:,3))-min(branchPointCluster(:,3)))*0.5 + (max(tmp(:,3))-min(tmp(:,3)))*0.5;      
    branchPointCluster(:,1) = branchPointCluster(:,1) - xdiff;
    branchPointCluster(:,2) = branchPointCluster(:,2) - ydiff;
    branchPointCluster(:,3) = branchPointCluster(:,3) - zdiff;   
end
function retTable = write2tablenew(lasFile)
    retTable = zeros(size(lasFile,1),3);
    retTable(:,1) = lasFile.x;
    retTable(:,2) = lasFile.y;
    retTable(:,3) = lasFile.z;
end
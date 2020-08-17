function VDAtt = getVoxelDimension(inputData, voxelSize)
    VDAtt.voxelSize = voxelSize;
    VDAtt.xDiv = round((inputData.XMax - inputData.XMin)/voxelSize,0);
    VDAtt.yDiv = round((inputData.YMax - inputData.YMin)/voxelSize,0);
    VDAtt.zDiv = round((inputData.ZMax - inputData.ZMin)/(voxelSize),0);
    VDAtt.xStep = (inputData.XMax - inputData.XMin)/VDAtt.xDiv;
    VDAtt.yStep = (inputData.YMax - inputData.YMin)/VDAtt.yDiv;
    VDAtt.zStep = (inputData.ZMax - inputData.ZMin)/VDAtt.zDiv;
    
    % get x,y, and z voxel divisions array
    [VDAtt.rX,VDAtt.rY,VDAtt.rZ] = getVoxelDivisions(inputData, VDAtt.xDiv, VDAtt.yDiv, VDAtt.zDiv);
end
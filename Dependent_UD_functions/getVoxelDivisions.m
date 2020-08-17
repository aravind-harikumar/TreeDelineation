function[rX,rY,rZ] = getVoxelDivisions(pcData, xDiv, yDiv, zDiv)
% Perform division space using: [(minX + size_of_one_division_along_X) * scale_vector_of_size_xDiv]
rX = pcData.XMin + (pcData.XMax - pcData.XMin)/(xDiv-1) * [0:xDiv];
rY = pcData.YMin + (pcData.YMax - pcData.YMin)/(yDiv-1) * [0:yDiv];
rZ = pcData.ZMin + (pcData.ZMax - pcData.ZMin)/(zDiv-1) * [0:zDiv];
end
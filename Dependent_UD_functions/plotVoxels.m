
function [F,V,C] = plotVoxels(Indices, M, ipc, printVoxel, voxelTransparency, colorArr)

        [F,V,C]=ind2patch(Indices,M,'v'); %Creating patch data for selection of low voxels
        
        % Scaling and cetering
        minV1  = min(V(:,1));  maxV1 = max(V(:,1));
        minV2  = min(V(:,2));  maxV2 = max(V(:,2));
        minV3  = min(V(:,3));  maxV3 = max(V(:,3));
        V(:,1) = (V(:,1) - (minV1))*(ipc.XMax - ipc.XMin)/(maxV1 - minV1) - (abs(ipc.XMin)); % origDat - shifting*scaling - centering
        V(:,2)=  (V(:,2) - (minV2))*(ipc.YMax - ipc.YMin)/(maxV2 - minV2) - (abs(ipc.YMin)); % origDat - shifting*scaling - centering
        V(:,3)=  (V(:,3) - (minV3))*(ipc.ZMax - ipc.ZMin)/(maxV3 - minV3) - (abs(ipc.ZMin)); % origDat - shifting*scaling - centering

        if(printVoxel)
            hs=patch('Faces',F,'Vertices',V,'EdgeColor','k', 'FaceVertexCData', repmat(colorArr, size(F,1), 1), 'FaceColor','flat','FaceAlpha',voxelTransparency);
        end
end
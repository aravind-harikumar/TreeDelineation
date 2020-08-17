function [vCellData, PCData] = formVoxelSpacefromTLSProcessing(PCData, VDim, printVoxel)

%OT = OcTree(PCData.lidarDataArray(:,1:3));
[vCellData, PCData] = getVoxelDensityAndNebrCellIndicesNew(PCData, VDim.rX, VDim.rY, VDim.rZ, VDim.xDiv, VDim.yDiv, VDim.zDiv);

hold on;
if(printVoxel)
    
    cnt=1;
    for jj = 1:1:1
        %uip = uipanel('Position',[0.75 0.75 0.25 0.25]);
        %subplot(1,2,cnt);
        Indices2 = find(vCellData.voxNeighArrDensity(:)>=jj);
        [M1,~,~] = meshgrid(VDim.rX(1:end-1),VDim.rY(1:end-1),VDim.rZ(1:end-1));
        [F,V,C]=ind2patch(Indices2,M1,'v'); %Creating patch data for selection of low voxels
        % Scaling and cetering
        minV1  = min(V(:,1));  maxV1 = max(V(:,1));
        minV2  = min(V(:,2));  maxV2 = max(V(:,2));
        minV3  = min(V(:,3));  maxV3 = max(V(:,3));
        V(:,1) = (V(:,1) - (minV1))*(PCData.XMax - PCData.XMin)/(maxV1 - minV1) - (abs(PCData.XMin)); % origDat - shifting*scaling - centering
        V(:,2)=  (V(:,2) - (minV2))*(PCData.YMax - PCData.YMin)/(maxV2 - minV2) - (abs(PCData.YMin)); % origDat - shifting*scaling - centering
        V(:,3)=  (V(:,3) - (minV3))*(PCData.ZMax - PCData.ZMin)/(maxV3 - minV3) - (abs(PCData.ZMin)); % origDat - shifting*scaling - centering
        voxelTransparency = 0.5;
        patch('Faces',F,'Vertices',V, 'EdgeColor',[0 1 0], 'CData', ones(size(C))*100,'FaceColor',[0 0.5 0],'FaceAlpha',voxelTransparency);
        axis equal;
        camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on;
        cnt = cnt + 1;
        axis([PCData.XMin PCData.XMax PCData.YMin PCData.YMax 0 ceil(PCData.ZMax)])
        setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera')
        zoom(1/4);
    end
end

end
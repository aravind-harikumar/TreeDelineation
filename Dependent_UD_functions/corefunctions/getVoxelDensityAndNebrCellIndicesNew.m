function [vCellData, PCData] = getVoxelDensityAndNebrCellIndicesNew(PCData, xedge, yedge, zedge, xDiv, yDiv, zDiv)

vCellData.voxIndxCellArr = cell(xDiv,yDiv,zDiv); % cell
vCellData.voxPCountArray= zeros(xDiv,yDiv,zDiv); %array
vCellData.voxNeighIndxCellArr = cell(xDiv,yDiv,zDiv); %cell
vCellData.voxelCenterCellArr = repmat({[0 0 0]},xDiv,yDiv,zDiv); %cell(xDiv,yDiv,zDiv); % cell
vCellData.voxelCenterCellXArr = zeros(xDiv,yDiv,zDiv);
vCellData.voxelCenterCellYArr = zeros(xDiv,yDiv,zDiv);
vCellData.voxelCenterCellZArr = zeros(xDiv,yDiv,zDiv);
vCellData.voxCellCornersCellArr = cell(xDiv,yDiv,zDiv); % cell
vCellData.voxContainedPointIndice = cell(xDiv,yDiv,zDiv); % cell
vCellData.voxNearCentroidPointIndex = zeros(xDiv,yDiv,zDiv); %array
vCellData.voxNeighArrDensity = repmat(-2,xDiv,yDiv,zDiv);

xlen = length(xedge)-1; xmid = zeros(xlen,1);
ylen = length(yedge)-1; ymid = zeros(ylen,1);
zlen = length(zedge)-1; zmid = zeros(zlen,1);
for i = 1:xlen
  xmid(i,1) = (xedge(i)+xedge(i+1))/2;
end
for i = 1:ylen
  ymid(i,1) = (yedge(i)+yedge(i+1))/2;
end
for i = 1:zlen
  zmid(i,1) = (zedge(i)+zedge(i+1))/2;
end

loc = zeros(size(PCData.lidarDataArray(:,1:3)));
[~,loc(:,1)] = histc(PCData.lidarDataArray(:,1),xedge);
[~,loc(:,2)] = histc(PCData.lidarDataArray(:,2),yedge);
[~,loc(:,3)] = histc(PCData.lidarDataArray(:,3),zedge);

[vCoord,IA,IC] = unique(loc,'rows');
PCData.lidarDataArray(:,9) = IC;
PCData.lidarDataArrayComplete(:,9) = IC;

for i = 1:1:size(vCoord,1)    
    %op = find(ismember(loc,vCoord(i,:),'rows'));
%    try
    selectedIndice = IC==IC(IA(i));  
    
    %vCellData.voxIndxCellArrnew{vCoord(i,1),vCoord(i,2),vCoord(i,3)} = op;
    vCellData.voxPCountArray(vCoord(i,1),vCoord(i,2),vCoord(i,3)) = sum(selectedIndice);
%     catch(Exception e)
%         sdfs = 0;
%     end
    voxelDataSubset = PCData.lidarDataArray(selectedIndice,:);
    vCellData.voxContainedPointIndice{vCoord(i,1),vCoord(i,2),vCoord(i,3)} = voxelDataSubset(:,9);
    
    xstep = xedge(vCoord(i,1)); ystep = yedge(vCoord(i,2)); zstep = zedge(vCoord(i,3));
    xindx = vCoord(i,1); yindx = vCoord(i,2); zindx = vCoord(i,3);    
    vCellData.voxelCenterCellArr{vCoord(i,1),vCoord(i,2),vCoord(i,3)} = [(xstep + xstep+abs(xstep-xedge(xindx+1)))/2 (ystep + ystep+abs(ystep-yedge(yindx+1)))/2 (zstep + zstep+abs(zstep-zedge(zindx+1)))/2 ];  %mean(voxelDataSubset,1);
    vCellData.voxelCenterCellXArr(vCoord(i,1),vCoord(i,2),vCoord(i,3)) = (xstep + xstep+abs(xstep-xedge(xindx+1)))/2;  %mean(voxelDataSubset,1);
    vCellData.voxelCenterCellYArr(vCoord(i,1),vCoord(i,2),vCoord(i,3)) = (ystep + ystep+abs(ystep-yedge(yindx+1)))/2;
    vCellData.voxelCenterCellZArr(vCoord(i,1),vCoord(i,2),vCoord(i,3)) = (zstep + zstep+abs(zstep-zedge(zindx+1)))/2;
   
    vCellData.voxCellCornersCellArr{vCoord(i,1),vCoord(i,2),vCoord(i,3)} = [xstep ystep zstep xstep+abs(xstep-xedge(xindx+1)) ystep+abs(ystep-yedge(yindx+1)) zstep+abs(zstep-zedge(zindx+1))];
    
    voxelIndx = sub2ind(size(vCellData.voxPCountArray),vCoord(i,1),vCoord(i,2),vCoord(i,3));
    vCellData.voxIndxCellArr{vCoord(i,1),vCoord(i,2),vCoord(i,3)} = voxelIndx;
    
    PCData.lidarDataArray(selectedIndice,9) = voxelIndx;
    PCData.lidarDataArray(selectedIndice,10:12) =  repmat(vCellData.voxelCenterCellArr{vCoord(i,1),vCoord(i,2),vCoord(i,3)}, sum(selectedIndice) ,1 );
    
    centroidPoint = mean(voxelDataSubset(:,1:3),1);
    [~,Inds] = pdist2(voxelDataSubset(:,1:3),centroidPoint(:,1:3),'euclidean','Smallest',1); % index of closest point to centrod
    if(~isempty(Inds))
        vCellData.voxNearCentroidPointIndex(vCoord(i,1),vCoord(i,2),vCoord(i,3)) = voxelDataSubset(Inds,7); % 
    end
    
    [Iadj, ea, ~] = neighbourND(voxelIndx, [xDiv yDiv zDiv], [0.1 0.1 0.1]);
    if( isequal( sum(ea == 0.1), 6 ))
        Iadj = Iadj(ea == 0.1);  
    else
        Iadj = [];
    end
    
    vCellData.voxNeighIndxCellArr{vCoord(i,1),vCoord(i,2),vCoord(i,3)} = Iadj;
    
end
voxIndWithPoints = find(vCellData.voxPCountArray>0);
% vCellData.voxArrDensity
% vCellData.voxNeighIndxCellArr
for i = 1:1:size(vCoord,1)  
     indx = sub2ind(size(vCellData.voxPCountArray),vCoord(i,1),vCoord(i,2),vCoord(i,3));
     if(ismember(indx,voxIndWithPoints))
        vCellData.voxNeighArrDensity(vCoord(i,1),vCoord(i,2),vCoord(i,3)) = sum(vCellData.voxPCountArray(vCellData.voxNeighIndxCellArr{vCoord(i,1),vCoord(i,2),vCoord(i,3)})>0);
     end
  %   [in, ff] =  knnsearch(vCellData.voxelCenterCellArr,vCellData.voxelCenterCellArr{44},'dist','cityblock','k',4);
         
end

end
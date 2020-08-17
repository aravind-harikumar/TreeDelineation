  function unAssVoxelIndxAndClosestAssignedINdx = temp(xDiv, yDiv, zDiv, voxelWithClassAssndIndx, searchIndices)
     [X,Y,Z] = meshgrid(linspace(1,xDiv,xDiv),linspace(1,yDiv,yDiv),linspace(1,zDiv,zDiv));
     points = [X(:) Y(:) Z(:)];
     voxelWithClassAssnd = points(voxelWithClassAssndIndx,:);     
     searchVoxels = points(searchIndices,:);
     [~,I] = pdist2(voxelWithClassAssnd,searchVoxels,'euclidean','Smallest',1);  
     unAssVoxelIndxAndClosestAssignedINdx = [searchIndices  voxelWithClassAssndIndx(I)];
end
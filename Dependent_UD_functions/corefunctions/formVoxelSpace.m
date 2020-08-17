function formVoxelSpace(Indices, projData, VDim, printVoxel, colorArr)
              
        % Form voxels
        [M,~,~]=meshgrid(VDim.rY(1:end-1),VDim.rX(1:end-1),VDim.rZ(1:end-1));     
        SubdomVoxelsByTreeArr = {};
        for treeCnt = 1:1:size(Indices,3)
            SubdomVoxelsByTree = zeros(size(M));
            for i = 1:1:size(Indices,4)
                SubdomVoxelsByTree(:,:,i) = flipud(Indices(:,:,treeCnt,i));
            end
            SubdomVoxelsByTreeArr{treeCnt} = SubdomVoxelsByTree;
        end
        
        F = 0; V =0; C =0;
        voxelTransparency = 0.5;  
        if(and(NC.ISPLOTON,printVoxel))
            f5 = figure;     
            set(f5, 'Position', [1320 60 600 420])
 
            A1 = SubdomVoxelsByTreeArr{1};
            Indices11 = find(A1(:)==200);
            Indices11 = [Indices11; 1; VDim.xDiv*VDim.yDiv*VDim.zDiv];
            [F,V,C] = plotVoxels(Indices11, M, projData, printVoxel, voxelTransparency, [0 0.4 0]);
            hold on;
            
            for i = 1:1:size(SubdomVoxelsByTreeArr,2) 
                A1 = SubdomVoxelsByTreeArr{i};
                Indices1 = [find(A1(:)==(20*i)); 1; VDim.xDiv*VDim.yDiv*VDim.zDiv];
                
                [F,V,C] = plotVoxels(Indices1, M, projData, printVoxel, voxelTransparency, colorArr(i,:));
                hold on; camproj perspective; rotate3d on; axis vis3d; axis on; grid on;
                view(135,60); box on; axis square;
            end
            axis([projData.XMin projData.XMax projData.YMin projData.YMax projData.ZMin projData.ZMax]);
        end
        
end

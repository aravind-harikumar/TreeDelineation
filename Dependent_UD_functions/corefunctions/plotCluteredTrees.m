function [getTreeTops, X_Extend, Y_Extend] = plotCluteredTrees(stData, pdData, subTreeVoxelsByFeature, indxArr, plotResult,colorArr)

    getTreeTops = []; X_Extend=[]; Y_Extend=[]; 
    
    % Adding new columns to stData.lidarDataArray to stored point class.
    stData.lidarDataArray = [stData.lidarDataArray zeros(size(stData.lidarDataArray,1),1)];

    % Plot the candidate segment (top left plot)
    if(and(NC.ISPLOTON,plotResult))
%         figure;
        subplot(2,4,[1 2]);
        plot3(stData.lidarDataArray(:,1), stData.lidarDataArray(:,2), stData.lidarDataArray(:,3),'.', 'MarkerSize', 10,'Color',[0 0.5 0]);
        xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
        camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45,45); box on;
        title( strcat('Candidate Segment Point Cloud'));
    end

    % Voxel cloud of Subdominant trees obtained in the projected space
    SubdomVoxelsByTreeArr = {};
    for treeCnt = 1:1:size(subTreeVoxelsByFeature,3)
        SubdomVoxelsByTree = zeros(size(subTreeVoxelsByFeature,1),size(subTreeVoxelsByFeature,2),size(subTreeVoxelsByFeature,4));
        for i = 1:1:size(subTreeVoxelsByFeature,4)     
            SubdomVoxelsByTree(:,:,i) = subTreeVoxelsByFeature(:,:,treeCnt,i);
        end
        SubdomVoxelsByTreeArr{treeCnt} = SubdomVoxelsByTree;
    end
    
    % Map-back the subdominant tree points to the original space
    TreeCountMax = size(subTreeVoxelsByFeature,3);
    subtreeData = {};
    AllAssignedPointArr = [];
    
    plotSlot = [3 4 7 8];  
    for iSubTreeCount = 1:1:TreeCountMax 
        subTreeVoxels = SubdomVoxelsByTreeArr{iSubTreeCount};
        subTreeVoxelIndices = find(subTreeVoxels==iSubTreeCount*20);
        
        pointIndicesInVoxelCell = ismember(pdData.lidarDataArray(:,9), subTreeVoxelIndices);        
        %Assign the indices to the 10th column of the lidarDataArray
        stData.lidarDataArray(pointIndicesInVoxelCell,8) = iSubTreeCount;
        %Store the assigned point indices in an array.
        AllAssignedPointArr = [AllAssignedPointArr; find(pointIndicesInVoxelCell==1)];
        
        VoxelCellPointData = ismember(stData.lidarDataArray(:,6), pdData.lidarDataArray(pointIndicesInVoxelCell,8));
        subtreeData{iSubTreeCount} = stData.lidarDataArray(VoxelCellPointData,1:3);
        %subtreeData = ExpandCloud(subtreeData, 1);  % Scaling factor Here is 1 --> actual size (for testing only)
    end
    % get dominant tree points
    [retMaxXYZ, ~] = findMaxHeightXY(stData.lidarDataArray(:,1:3)); 
    in = incircle(stData.lidarDataArray(:,1:2),retMaxXYZ(1,1:2), 1);
    stData.lidarDataArray(in,8) = -1;
    
    opd = stData.lidarDataArray(:,8)~=0;
    AssignedVoxelCellPointData = stData.lidarDataArray(opd,:);
    opdd = stData.lidarDataArray(:,8)==0;
    NonAssignedVoxelCellPointData = stData.lidarDataArray(opdd,:);
    
    [dist,Indx] = pdist2(AssignedVoxelCellPointData(:,1:3),NonAssignedVoxelCellPointData(:,1:3),'euclidean', 'Smallest' ,1);
    NonAssignedVoxelCellPointData(:,8) = AssignedVoxelCellPointData(Indx,8);
    stData.lidarDataArray = [NonAssignedVoxelCellPointData; AssignedVoxelCellPointData];
    
    % plot 3D data of subdominat trees
    sTreeIndicesArr = unique(stData.lidarDataArray(:,8));
    % restrict to displaying only first 4 trees
    subtreeCount = size(sTreeIndicesArr,1)-1;
    if(subtreeCount >4)
        subtreeCount = 4;    end
    for iSubTreeCount = 1:1:subtreeCount
        indx = stData.lidarDataArray(:,8) == sTreeIndicesArr(iSubTreeCount+1); % first indx is dom tree data
        subtreeData = stData.lidarDataArray(indx,:);
        if(~isempty(subtreeData))
            if(and(NC.ISPLOTON,plotResult))
                subplot(2,4,plotSlot(iSubTreeCount));
                plot3(subtreeData(:,1), subtreeData(:,2), subtreeData(:,3), '.','MarkerSize', 10,'Color',colorArr(iSubTreeCount,:));
                hold on;
                xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
                camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45,45); box on;
                axis([min(subtreeData(:,1))-1 max(subtreeData(:,1))+1 min(subtreeData(:,2))-1 max(subtreeData(:,2))+1 0 ceil(stData.treeHeight)])
                title( strcat({'Subdominant Tree '}, num2str(iSubTreeCount)));
            end
        end
        
        if(sum(indx)>0)
            getTreeTops = [getTreeTops; [mean(stData.origData(VoxelCellPointData,1)) mean(stData.origData(VoxelCellPointData,2))] ];
            X_Extend = [X_Extend; [max(stData.origData(VoxelCellPointData,1)) min(stData.origData(VoxelCellPointData,1))]];
            Y_Extend = [Y_Extend; [max(stData.origData(VoxelCellPointData,2)) min(stData.origData(VoxelCellPointData,2))]];
        else
            getTreeTops = [0 0]; X_Extend = [0 0];  Y_Extend = [0 0];
        end
    end
    hold on;
    % plot 3D data of the dominant tree
    indx = stData.lidarDataArray(:,8) == -1;
    subtreeData = stData.lidarDataArray(indx,:);
    if(and(NC.ISPLOTON,plotResult))
        subplot(2,4,[5 6]);
        plot3(subtreeData(:,1), subtreeData(:,2), subtreeData(:,3), '.','MarkerSize', 10,'Color',[0 0.5 0]);

        xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
        camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45,45); box on;
        axis([min(subtreeData(:,1))-1 max(subtreeData(:,1))+1 min(subtreeData(:,2))-1 max(subtreeData(:,2))+1 0 ceil(stData.treeHeight)])
        title( strcat({'Subdominant Tree '}, num2str(iSubTreeCount)));
    end
    
%     % Plot the first 4 subdominat trees in the original space 
%     AllSubTreeIndex = []; getTreeTops = []; X_Extend=[]; Y_Extend=[]; 
%     for iSubTreeCount = 1:1:4 % Set TreeCountMax = 4
%         
%         if(and(NC.ISPLOTON,plotResult))
%             subplot(2,4,plotSlot(iSubTreeCount));    
%         end
%         
%         if(size(subtreeData,1)>0)        
%             in = inhull(stData.lidarDataArray(:,1:3) ,subtreeData(:,1:3), convhulln(subtreeData(:,1:3)));
%             AllSubTreeIndex = [AllSubTreeIndex; stData.lidarDataArray(in,6)];
% 
%             if(sum(VoxelCellPointData)>0)
%                 getTreeTops = [getTreeTops; [mean(stData.origData(VoxelCellPointData,1)) mean(stData.origData(VoxelCellPointData,2))] ];
%                 X_Extend = [X_Extend; [max(stData.origData(VoxelCellPointData,1)) min(stData.origData(VoxelCellPointData,1))]];
%                 Y_Extend = [Y_Extend; [max(stData.origData(VoxelCellPointData,2)) min(stData.origData(VoxelCellPointData,2))]];
% 
%                 if(and(NC.ISPLOTON,plotResult))
%                     %plot3(treeData(:,1), treeData(:,2), treeData(:,3),'.', 'MarkerSize', 10,'Color',[0 0.5 0]);
%                     %hold on;
%                     plot3(stData.lidarDataArray(in,1), stData.lidarDataArray(in,2), stData.lidarDataArray(in,3), '.','MarkerSize', 10,'Color',[0 0.5 0]);
% 
%                     xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
%                     camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45,45); box on;
%                     axis([min(subtreeData(:,1))-1 max(subtreeData(:,1))+1 min(subtreeData(:,2))-1 max(subtreeData(:,2))+1 0 ceil(stData.treeHeight)])
%                     title( strcat({'Subdominant Tree '}, num2str(iSubTreeCount)));
%                 end
%             end
%         end
%         
%     end
%     
%     % Plot the dominant tree
%     optree = ~ismember(stData.lidarDataArray(:,6), unique(AllSubTreeIndex));
%     
%     accumIndxArr = [];
%     subTreeVoxels = SubdomVoxelsByTreeArr{1};
%     subTreeVoxelIndices = find(subTreeVoxels==200);    
%     for j=1:1:size(subTreeVoxelIndices,1);
%         accumIndxArr = [accumIndxArr; indxArr{subTreeVoxelIndices(j)}];
%     end
%     pointIndicesInVoxelCell = ismember(stData.index, accumIndxArr);
%     if(and(NC.ISPLOTON,plotResult))
%         subplot(2,4,[5 6]);
%         plot3(stData.lidarDataArray(optree,1), stData.lidarDataArray(optree,2), stData.lidarDataArray(optree,3), '.', 'MarkerSize', 10, 'Color', [0 0.5 0]);
%         xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
%         camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45,45); box on;
%         title('Dominant Tree');
%     end
    
end


% 
%  op = ismember(pdData.lidarDataArray(:,9), subTreeVoxelIndices);
%         opp = ismember(stData.lidarDataArray(:,6),pdData.lidarDataArray(op,8));       
%         treeData = stData.lidarDataArray(opp,1:3);
%         
%         treeData(:,3) = -treeData(:,3) + max(treeData(:,3));
%         
%         tmpForHull = ExpandCloud(treeData, 1.2);
%         
%         in = inhull(stData.lidarDataArray(:,1:3) ,tmpForHull(:,1:3), convhulln(tmpForHull(:,1:3)));
% %         
%          xShift =  mean(treeData(:,1)); 
%          yShift =  mean(treeData(:,2)); 
%          maxZVal = max(treeData(:,3));
% 
%          isParabloid = true; pltParaboloid = true; % change if another geometric shape is used
%          coeffArr = fitParabloid(tmpForHull);
%          plotShape(coeffArr, xShift, yShift, maxZVal, isParabloid, pltParaboloid);
%        
%          hold on;
%          AllSubTreeIndex = [AllSubTreeIndex; stData.lidarDataArray(in,6)];
%        
%          aa = [];
%          slabbedCoordinates1 = -stData.lidarDataArray(:,3) + max(stData.lidarDataArray(:,3));
%          slabbedCoordinates1 = slabbedCoordinates1(and(abs(slabbedCoordinates1(:,1))<2.5, abs(slabbedCoordinates1(:,2))<2.5),:);
%          for pt =1:1:size(slabbedCoordinates1,1);                    
%             isIncluded = isPointIncluded(slabbedCoordinates1(pt,1:3),coeffArr(1),coeffArr(2),coeffArr(3),1,'elparaboloid');
%             if(~isIncluded)
%                 % plot3(slabbedCoordinates1(pt,1), slabbedCoordinates1(pt,2), slabbedCoordinates1(pt,3),'.', 'MarkerSize' ,40, 'Color',[0 0 0.5]);
%                 aa = [aa slabbedCoordinates1(pt,4)];
%             end
%          end 
%         
%         if(sum(opp)>0)
%             getTreeTops = [getTreeTops; [mean(stData.origData(opp,1)) mean(stData.origData(opp,2))] ];
%             X_Extend = [X_Extend; [max(stData.origData(opp,1)) min(stData.origData(opp,1))]];
%             Y_Extend = [Y_Extend; [max(stData.origData(opp,2)) min(stData.origData(opp,2))]];
%             
%             if(and(NC.ISPLOTON,plotResult))
%                 %plot3(treeData(:,1), treeData(:,2), treeData(:,3),'.', 'MarkerSize', 10,'Color',[0 0.5 0]);
%                 %hold on;
%                 plot3(stData.lidarDataArray(in,1), stData.lidarDataArray(in,2), stData.lidarDataArray(in,3), '.','MarkerSize', 10,'Color',[0 0.5 0]);
%         
%                 xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
%                 camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45,45); box on;
%                 axis([min(treeData(:,1)) max(treeData(:,1)) min(treeData(:,2)) max(treeData(:,2)) 0 ceil(stData.treeHeight)])
%                 title( strcat({'Subdominant Tree '}, num2str(iSubTreeCount)));
%             end 
%             
%         end
%     end
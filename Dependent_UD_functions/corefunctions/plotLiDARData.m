function plotLiDARData(stData, plotStem) 

    if(and(NC.ISPLOTON,plotStem))        
        f1=figure(1);
        set(f1, 'Position', [10 620 950 500])
        subplot(1,2,1);
        hold on;
        plotXScale = -10:1:10; plotYScale = -10:1:10;
        plot3([0;0], [0;0], [0;stData.retMaxXYZ(3)],'-*', 'MarkerSize', 20);
    end
   
    maxXYZ = stData.retMaxXYZ;
    retMaxPointIndex = stData.retMaxPointIndexNew;
        
    org = stData.origData;    
    aa = org(:,1:3);
    
    PCAplotOn = false;
    REFpoint  = getPCAVectors(stData,PCAplotOn);
    opOffseted = [REFpoint(1:2,:) [10;10]];    
    bb = repmat(opOffseted(1,:),size(org,1),1);

    [~, cc]= min(abs(sqrt(sum((aa.^2-bb.^2), 2))));   
    indexExtreme = cc;   
    indexTreeTop = stData.retMaxXYZIndex;
    retTreeTop = [stData.lidarDataArray(indexTreeTop,1) stData.lidarDataArray(indexTreeTop,2) stData.lidarDataArray(indexTreeTop,3)];
      
    ppa = [retTreeTop(1) retTreeTop(2) 0]; 
    ppb =  opOffseted(1,:); %[stData.lidarDataArray(indexExtreme,1) stData.lidarDataArray(indexExtreme,2) stData.lidarDataArray(indexExtreme,3)];
    ppc = retTreeTop;    
    planePointsAll = [ppa;ppb;ppc;ppa];
    
    if(NC.ISPLOTON)
        plot3(planePointsAll(3,1), planePointsAll(3,2) ,planePointsAll(3,3), '.', 'MarkerSize', 40,'Color', [1 0 0]);
        %text(ppa(1),ppa(2),ppa(3),'Projection of tree top on XY plane','HorizontalAlignment','right');
        %text(ppb(1),ppb(2),ppb(3),'reference point','HorizontalAlignment','right');
        %text(ppc(1),ppc(2),ppc(3),'Tree top', 'HorizontalAlignment','right');
        camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(-45, 15); 
        axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
        %plot3(stData.lidarDataArray(:,1), stData.lidarDataArray(:,2), stData.lidarDataArray(:,3),'.','MarkerSize', 10,'Color',[0 0.5 0]);
        
        %stData.lidarDataArray = [stData.lidarDataArray; 9 9 1 0 0 0 0]; 
        
        abc = [stData.lidarDataArray(:,1) stData.lidarDataArray(:,2)];
        
        op =  pdist2(abc,[0 0],'euclidean');

        uu = scatter3(stData.lidarDataArray(:,1), stData.lidarDataArray(:,2), stData.lidarDataArray(:,3),200, op,'.');
        hold on; grid on;
        
        % plot 3D plane
        %plotPlaneFrom3Points(ppa, ppb, ppc, true);

        hold off;
        ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 30);
        xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 30);
        zlabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 30);
        %title('LiDAR Data');
        set(gca,'fontsize',30)
        %printVoxel = false; printAll = false; voxelTransparency = 1;
        %hs = formVoxelSpace(stData.lidarDataArray(:,1:3), xDiv, yDiv, zDiv, printAll, printVoxel, voxelTransparency);
        colormap(summer);
        sdatatemp = stData.lidarDataArray;
        axis([-7 7 -7 7 0 25]);
        
        %grid minor
        %axis([min(sdatatemp(:,1)) max(sdatatemp(:,1)) min(sdatatemp(:,2))...
        %max(sdatatemp(:,2)) min(sdatatemp(:,3)) max(sdatatemp(:,3)) ]);
    end
 
end

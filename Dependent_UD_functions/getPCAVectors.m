function [OfsettedPC1Intersects, OfsettedPC2Intersects, OfsettedResIntersects,DiffDist, EV_ratio_thresh, EVals] = getPCAVectors(stData,plotOn)

        indexTreeTop = stData.retMaxXYZIndex;
        retTreeTop = [stData.lidarDataArray(indexTreeTop,1) stData.lidarDataArray(indexTreeTop,2) stData.lidarDataArray(indexTreeTop,3)];      
        ppc = [retTreeTop(1) retTreeTop(2)];

        %stData.lidarDataArray = 0;
        %stData.lidarDataArray = [0 0; 10 0; -10 0; -5 0; 5 0; 0 -18;];
        
        [COEFF, SCORE, LATENT, ~] = pca(stData.lidarDataArray(:,1:2));
        EVals = LATENT(1:2);
        
       
        EV_ratio_thresh = 0.5;
        if(and(LATENT(1)>0, LATENT(2)>0))
            EV_ratio_thresh = LATENT(2)/LATENT(1); %  % need to vary based on desnity and sub-tree height is the revised version.
        end
%         EV_ratio_thresh = 0.67;
        
        centeredData = SCORE*COEFF';
        PCVector = COEFF.*repmat(LATENT',2,1);
        
        % resultant direction
        Resultant = getVectorResultant(PCVector(:,1)',PCVector(:,2)');
        
        if(plotOn)
            figure;
            plot(centeredData(:,1), centeredData(:,2), '.','MarkerSize', 10,'Color',[0 0.5 0]);
            plot(ppc(:,1), ppc(:,2), '.', 'MarkerSize', 40,'Color', [1 0 0]);
            text(ppc(1),ppc(2),'Tree top', 'HorizontalAlignment','right');
            plotv(PCVector,'-r');
            hold on;
            plotv(-PCVector,'-r');
            text(PCVector(1,1),PCVector(2,1),'PC1', 'HorizontalAlignment','right');
            text(PCVector(1,2),PCVector(2,2),'PC2', 'HorizontalAlignment','right');
        end
        % get boundary point in reference direction
        k = boundary(centeredData(:,1),centeredData(:,2),0.5);
        boundaryPoints = [centeredData(k,1) centeredData(k,2)];
        pca1Direction = PCVector(:,2)';
        pca1Direction = [pca1Direction*100; pca1Direction*-100];
        
        pca2Direction = PCVector(:,1)';
        pca2Direction = [pca2Direction*100; pca2Direction*-100];
        
        Resultant = [Resultant*100; Resultant*-100];
        
        % intersection of pc1 with the XY projected cloud data boundary
        [X_pc1intersect,Y_pc1intersect,~,~] = intersections(pca1Direction(:,1),pca1Direction(:,2),boundaryPoints(:,1),boundaryPoints(:,2),'robust');
         
        % intersection of pc1 with the XY projected cloud data boundary
        [X_pc2intersect,Y_pc2intersect,~,~] = intersections(pca2Direction(:,1),pca2Direction(:,2),boundaryPoints(:,1),boundaryPoints(:,2),'robust');
        
        % intersection of resultant with the XY projected cloud data boundary
        [X_ResultantIntersect,Y_ResultantIntersect,~,~] = intersections(Resultant(:,1),Resultant(:,2),boundaryPoints(:,1),boundaryPoints(:,2),'robust');
        
        offsetXY =  centeredData(1,1:2)-stData.lidarDataArray(1,1:2); %first point
        
        % offsetted for returning pc1 value
        XPC1YPC1 = [X_pc1intersect Y_pc1intersect];
        OfsettedPC1Intersects =  XPC1YPC1 - repmat(offsetXY,size(XPC1YPC1,1),1);
        
        % offsetted for returning pc2 value
        XPC2YPC2 = [X_pc2intersect Y_pc2intersect];
        OfsettedPC2Intersects =  XPC2YPC2 - repmat(offsetXY,size(XPC2YPC2,1),1);
        
        % offsetted for returning value
        XResYRes = [X_ResultantIntersect Y_ResultantIntersect];
        OfsettedResIntersects =  XResYRes - repmat(offsetXY,size(XResYRes,1),1);
        
        
        % distance D1 between PC1-Boundary intersect point pairs
        distBtwPC1pairs = pdist(XPC1YPC1,'euclidean');
        % distance D2 between PC2-Boundary intersect point pairs
        distBtwPC2pairs = pdist(XPC2YPC2,'euclidean');
        % distance REDDIST between RES-Boundary intersect point pairs
        distBtwRESpairs = pdist(XResYRes,'euclidean');
        % ratio of (D1-D2)/D1;
        %DiffDist = abs(distBtwPC1pairs-distBtwRESpairs)/distBtwPC1pairs;
        DiffDist = min(distBtwPC1pairs,distBtwPC2pairs)/ max(distBtwPC1pairs,distBtwPC2pairs);
        
        if(plotOn)
            hold on;
            plot(boundaryPoints(:,1), boundaryPoints(:,2));
            plot(XPC1YPC1(:,1),XPC1YPC1(:,2),'.','MarkerSize', 20, 'Color','r');
            plot(XPC2YPC2(:,1),XPC2YPC2(:,2),'.','MarkerSize', 20, 'Color','b');
            plot(XResYRes(:,1),XResYRes(:,2),'.','MarkerSize', 20, 'Color','b');
            axis equal;
            hold on;
            plot(stData.lidarDataArray(:,1)+offsetXY(1), stData.lidarDataArray(:,2)+offsetXY(2), '.','MarkerSize', 3,'Color',[0.5 0.5 0.5]);
        end
        %plot(opOffseted(:,1),opOffseted(:,2),'*','MarkerSize', 10, 'Color','r');
        
%         axis([-10  10 -10 10]);
%         hold off; axis equal;
%         xlabel('Y Axis'); ylabel('X Axis');
%         title('LiDAR data of Tree');
            
end


% function [opOffseted, EV_ratio_thresh] = getPCAVectors(stData)
% 
%         indexTreeTop = stData.retMaxXYZIndex;
%         retTreeTop = [stData.lidarDataArray(indexTreeTop,1) stData.lidarDataArray(indexTreeTop,2) stData.lidarDataArray(indexTreeTop,3)];      
%         ppc = [retTreeTop(1) retTreeTop(2)];
% 
%         %stData.lidarDataArray = 0;
%         %stData.lidarDataArray = [0 0; 10 0; -10 0; -5 0; 5 0; 0 -18;];
%         
%         [COEFF, SCORE, LATENT, ~] = pca(stData.lidarDataArray(:,1:2));
%         
%         EV_ratio_thresh = 0.5;
%         if(and(LATENT(1)>0, LATENT(2)>0))
%             EV_ratio_thresh = LATENT(2)/LATENT(1);
%         end
%         
%         centeredData = SCORE*COEFF';
%         PCVector = COEFF.*repmat(LATENT',2,1);
%  %       figure(3);
%         
% %         plot(centeredData(:,1), centeredData(:,2), '.','MarkerSize', 10,'Color',[0 0.5 0]);
% %         plot(ppc(:,1), ppc(:,2), '.', 'MarkerSize', 40,'Color', [1 0 0]);
% %         text(ppc(1),ppc(2),'Tree top', 'HorizontalAlignment','right');
% %         plotv(PCVector,'-r')
% %         plotv(-PCVector,'-r')
% %         text(PCVector(1,1),PCVector(2,1),'PC1', 'HorizontalAlignment','right');
% %         text(PCVector(1,2),PCVector(2,2),'PC2', 'HorizontalAlignment','right');
%         
%         % get boundary point in reference direction
%         k = boundary(centeredData(:,1),centeredData(:,2),0.5);        
%         boundaryPoints = [centeredData(k,1) centeredData(k,2)];        
%         ab = PCVector(:,2)';
%         cd = PCVector(:,1)';
%         ab = [ab; ab*100; ab*-100];
%         cd = [cd; cd*100; cd*-100];
%         [x0,y0,~,~] = intersections(ab(:,1),ab(:,2),boundaryPoints(:,1),boundaryPoints(:,2),'robust');
%         [x00,y00,~,~] = intersections(cd(:,1),cd(:,2),boundaryPoints(:,1),boundaryPoints(:,2),'robust');
%         % offsetted for retuening value
%         XY = [x0 y0];
%         yyy = [x00 y00];   
%         
%         offsetXY =  centeredData(1,1:2)-  stData.lidarDataArray(1,1:2);
%         opOffseted =  XY - repmat(offsetXY,size(XY,1),1);
%         opOffsetedYY = yyy - repmat(offsetXY,size(yyy,1),1);
%         treetop = centeredData(indexTreeTop,1:2);
%         
%         ss = pdist2(treetop,[opOffseted(:,1:2);opOffsetedYY(:,1:2)],'Euclidean');
%         ss = sort(ss,'ascend');        
%         ssr = [ss(1)/ss(3) ss(1)/ss(4) ss(2)/ss(3) ss(2)/ss(4)];
%         [EV_ratio_thresh,indx] = min(ssr);
%         
% %         hold on;
% %         plot(treetop(:,1),treetop(:,2),'*','MarkerSize', 10, 'Color','k');
% %         hold on;
% %         plot(boundaryPoints(:,1), boundaryPoints(:,2));
% %         plot(XY(:,1),XY(:,2),'*','MarkerSize', 10, 'Color','r');
% %         axis equal;
% %         hold on;
%         %plot(stData.lidarDataArray(:,1), stData.lidarDataArray(:,2), '.','MarkerSize', 10,'Color',[0.5 0.5 0.5]);
%         %plot(opOffseted(:,1),opOffseted(:,2),'*','MarkerSize', 10, 'Color','r');
%         
% %         axis([-10  10 -10 10]);
% %         hold off; axis equal;
% %         xlabel('Y Axis'); ylabel('X Axis');
% %         title('LiDAR data of Tree');
% 
% end
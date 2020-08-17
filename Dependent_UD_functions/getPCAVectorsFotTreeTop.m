function opOffseted = getPCAVectorsFotTreeTop(stDatas)
       stData.lidarDataArray = stDatas;
       % indexTreeTop = stData.retMaxXYZIndex;
       % retTreeTop = [stData.lidarDataArray(indexTreeTop,1) stData.lidarDataArray(indexTreeTop,2) stData.lidarDataArray(indexTreeTop,3)];      
       % ppc = [retTreeTop(1) retTreeTop(2)];

       
        %plot(stData.lidarDataArray(:,1), stData.lidarDataArray(:,2), '.','MarkerSize', 10,'Color',[0.5 0.5 0.5]);
       
        [COEFF, SCORE, LATENT, ~] = pca(stData.lidarDataArray(:,1:2));
        
        centeredData = SCORE*COEFF';
        PCVector = COEFF.*repmat(LATENT',2,1);
        %figure(3); 
        hold on;
        plot(centeredData(:,1), centeredData(:,2), '.','MarkerSize', 10,'Color',[0 0.5 0]);
        %plot(ppc(:,1), ppc(:,2), '.', 'MarkerSize', 40,'Color', [1 0 0]);
        %text(ppc(1),ppc(2),'Tree top', 'HorizontalAlignment','right');
        
        plotv(PCVector,'-r')
        plotv(-PCVector,'-r')
        text(PCVector(1,1),PCVector(2,1),'PC1', 'HorizontalAlignment','right');
        text(PCVector(1,2),PCVector(2,2),'PC2', 'HorizontalAlignment','right');
        
  
        ab = PCVector(:,2)';        
        cd = PCVector(:,1)';        
        ab = [ab; ab*100; ab*-100];
        cd = [cd; cd*100; cd*-100];        
        [x0,y0,~,~] = intersections(ab(:,1),ab(:,2),cd(:,1),cd(:,2),'robust');
        XY = [x0 y0]; 
        offsetXY =  centeredData(1,1:2) - stData.lidarDataArray(1,1:2);
        opOffseted =  XY - repmat(offsetXY,size(XY,1),1);
        plot(opOffseted(:,1),opOffseted(:,2),'*','MarkerSize', 10, 'Color','r');        

        [csm, cdm] =  getSModel([centeredData stData.lidarDataArray(:,3)], 50, 50);
        
        
end
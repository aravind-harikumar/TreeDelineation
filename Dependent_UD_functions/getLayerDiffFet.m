function diffFetArr = getLayerDiffFet(slabbedCoordinates, numLayers, xDiv, yDiv, filterOrder)
    
    diffFetArr ={};
    zStep= max(slabbedCoordinates(:,3))/numLayers;    
    csm ={};
    % Generate layerwise DSMs (lowest layer first)
    for iLayerCount = numLayers:-1:1     
        cond1 = (slabbedCoordinates(:,3)> zStep*(iLayerCount-1));
        cond2 = (slabbedCoordinates(:,3)<= zStep*(iLayerCount));        
        slbCordByLayer = slabbedCoordinates(and(cond1,cond2),:);
        [csm{(numLayers+1)-iLayerCount}, ~] = getCSMCDM(slbCordByLayer,xDiv, yDiv, filterOrder, false, false);  
    end
    % calcualte differences of DSMs layers
    for jCount = 1:1:numLayers-1
        diffFetArr{jCount} = csm{jCount} - csm{jCount+1};
    end
end
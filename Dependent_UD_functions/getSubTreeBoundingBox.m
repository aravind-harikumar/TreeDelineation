function [crownBoundaryImage, colorArr] = getSubTreeBoundingBox(segmentedcsmcdm, csmNormalized,thresholdVal ,plotOn)
    
    I = mat2gray(csmNormalized);

    % I_ORIG, groundClearamnce, peakparam, envSmooth, plotYN;  
    [~, ~, crownBoundaryImage, colorArr, ~, ~, ~, ~] = ...
        subTreeTopDetectionInCSM(segmentedcsmcdm,I, thresholdVal, 0.0, 0.0, 1000, plotOn, 'ellipse'); % ellipse or rectangle
    
%     if(and(plotOn,NC.ISPLOTON))
%         figure;
%         imagesc(I);
%         hold on;
%        % imagesc(flip(crownBoundaryImage)); alpha(0.5);
%     end
 
end
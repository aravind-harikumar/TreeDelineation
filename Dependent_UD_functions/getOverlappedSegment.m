function overlapSegmArr = getOverlappedSegment(refImageOrg, glcmSliceImage, domTreeThrshLevel, levl, VOX_DIVS)        
    overlapSegmArr = zeros(VOX_DIVS.xDiv, VOX_DIVS.yDiv,size(refImageOrg,2));

    segCountIn = 0;
    for iRefSeg = 1:1:size(refImageOrg,2)  
        overlapSegm = zeros(size(glcmSliceImage));
        se = strel('disk',2,4);
        rimg = imresize(refImageOrg{iRefSeg},1);
        rimg(rimg>0) = 1;
        refImage = imdilate(rimg,se);
        % loop over glcm Slice segments
        
        imgClassArr = unique(glcmSliceImage);
        for iComp=1:1:length(unique(glcmSliceImage))           

            tempSlice = zeros(size(glcmSliceImage));
            tempSlice(glcmSliceImage==imgClassArr(iComp)) = 1;
            tempImg = double(tempSlice);
            
            % edges are cut off to avoid over segmentation problem.
            tempImg(1:size(refImageOrg,1),1:2)=0;
            tempImg(1:size(refImageOrg,1),end-2:end)=0;
            tempImg(end-2:end,1:size(refImageOrg,2))=0;
            
            % obtain overlap            
            indicesTemp = find(tempImg>0);
            indicesRef = find(refImage>0);
            
            overlappedMembers = ismember(indicesTemp, indicesRef);           
            nonoverlappedMembers = ~overlappedMembers;
            
            cond2 = true;% sum(nonoverlappedMembers)<(sum(refImage(:))*2.9);
            cond3 = true; %sum(nonoverlappedMembers)>0;
            cond1 = (sum(overlappedMembers(:))/size(overlappedMembers,1))>0.50;
            
            if(and(and(cond1,cond2),cond3))
                tempImg(indicesTemp(nonoverlappedMembers)) = 0;
                overlapSegm = overlapSegm + tempImg;
                segCountIn = segCountIn+1;
            end            
            
        end
        
        overlapSegm(overlapSegm==1) = 20*iRefSeg;
        overlapSegm = imdilate(overlapSegm,se);
        
        if(levl <=  round(domTreeThrshLevel,0) )
            %overlapSegm(overlapSegm==0) = 200;
            overlapSegm(overlapSegm>=0) = 200;
        else
           % levelOnFlag = false;
        end
        
       overlapSegmOld = overlapSegm;
        
        overlapSegm = imresize(overlapSegm, [VOX_DIVS.xDiv, VOX_DIVS.yDiv], 'nearest');
        overlapSegmArr(:,:,iRefSeg) = flipud(overlapSegm);
        
    end
    
end
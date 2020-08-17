function [csm, cdm] =  getSModel(slabCord, xDiv, yDiv)

    MMS = getMaxMins(slabCord);

    gXStep = abs(MMS.maxX - MMS.minX)/(xDiv);
    gYStep = abs(MMS.maxY - MMS.minY)/(yDiv);   

    halfXGridSize  = gXStep/2;       
    halfYGridSize  = gXStep/2; 

    xDimGrid = size(MMS.minX+halfXGridSize : gXStep : MMS.maxX-halfXGridSize,2);
    yDimGrid = size(MMS.minY+halfYGridSize : gYStep : MMS.maxY-halfYGridSize,2);        

    csm = zeros(xDimGrid,yDimGrid);
    cdm = zeros(xDimGrid,yDimGrid);

    cnt1 = 1; 

    for i = MMS.minX+halfXGridSize:gXStep:MMS.maxX-halfXGridSize
        cnt2 = 1;
        for j = MMS.minY+halfYGridSize:gYStep:MMS.maxY-halfYGridSize               
            cond1 = and(slabCord(:,1)>i-halfXGridSize, slabCord(:,1)<=i+halfXGridSize);
            cond2 = and(slabCord(:,2)>j-halfYGridSize, slabCord(:,2)<=j+halfYGridSize);
            indGridPoints = find(and(cond1,cond2));

            if(size(indGridPoints,1)>0)  
                csm(cnt1,cnt2) = mean(slabCord(indGridPoints,3));   
                cdm(cnt1,cnt2) = length(slabCord(indGridPoints,3));  
            end

            cnt2 = cnt2+1;
        end
        cnt1 = cnt1+1;
   end

   csm = ((csm-min(csm(:)))/( max(csm(:)) - min(csm(:)) ))*255;
   
   
   cdm = ((cdm-min(cdm(:)))/( max(cdm(:)) - min(cdm(:)) ))*255;
   
   imagesc(csm.*cdm);
   hold on;
   aff = mat2gray(csm.*cdm);
   temp = zeros(size(aff));   
   sdf = find(aff(:)>=(max(aff(:))-0.5));
   temp(sdf) = 10;   
   imshow(temp,[0 1]);
   alpha(0.5)
   
end
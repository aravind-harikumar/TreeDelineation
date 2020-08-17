function seeds = subTreePeaks(img, smoothThresh)
    imgNRows =size(img,1);
    img(img>0) = 1; img((imgNRows/2)+5:imgNRows,:) = 1;
    edgeImg = edge(img); 
    %imagesc(edgeImg);    
    [col,row] = ind2sub(size(edgeImg),find(edgeImg(:)>0));
    col = smooth(row,col,smoothThresh,'rloess');
    [loc,~] = peakfinder(imgNRows-col);
    seeds = [row(loc) col(loc)+((size(img,1)-col(loc))/2)];
    %hold on; plot(row(loc),col(loc)+((size(img,1)-col(loc))/2),'*','Color',[1 0 0 ]);
end
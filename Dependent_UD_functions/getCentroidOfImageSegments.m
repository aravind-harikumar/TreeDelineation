function centrods = getCentroidOfImageSegments(refImage)
    se = strel('disk',2,4);
    refImage = imdilate(refImage,se);
    BW = imbinarize(imresize(refImage,1));
    [B,~] = bwboundaries(BW,'noholes');
    centrods = [];
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
        [geom, ~, ~] = polygeom( boundary(:,2), boundary(:,1) );   
        centrods = [centrods; [geom(2) geom(3)]];
    end
end
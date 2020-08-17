function L = watershedNew(I, seedsIndx)
    %I =1-I;
    I = imadjust(I,[0 0.9],[]);
    
    I = imresize(I,2);
    %I = imgaussfilt(I,2);    
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    imshow(gradmag,[])
    %title('Gradient magnitude (gradmag)')    
    se = strel('disk', 3,4);
    Io = imopen(I, se);
    imagesc(Io),
    %title('Opening (Io)')
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    imagesc(Iobr)
    %title('Opening-by-reconstruction (Iobr)')
    Ioc = imclose(Io, se);
    imagesc(Ioc)
    %title('Opening-closing (Ioc)')
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    imagesc(Iobrcbr)
    %title('Opening-closing by reconstruction (Iobrcbr)')
    %IobrcbrCopy =Iobrcbr;
    %IobrcbrCopy(IobrcbrCopy<0.6) = 0;
    %Iobrcbr = imdilate(Iobrcbr, se);
    %IobrcbrCopy(IobrcbrCopy>0.8) =1;
    
    seedsIndx = sub2ind(size(I), seedsIndx(:,2), seedsIndx(:,1));
    seedsIndx = uint16(seedsIndx);
    seeds = zeros(size(Iobrcbr));
    seeds(seedsIndx) = 1;
    seeds = imdilate(seeds, se);
    
    fgm = imregionalmax(Iobrcbr,4);
    imagesc(fgm)
    %title('Regional maxima of opening-closing by reconstruction (fgm)')
    I2 = I;
    I2(fgm) = 1;
    imagesc(I2)
    %title('Regional maxima superimposed on original image (I2)')
    se2 = strel(ones(3,3));
    fgm2 = imclose(fgm, se2);
    fgm3 = imerode(fgm2, se2);
    fgm4 = bwareaopen(fgm3, 1);
    I3 = I;
    I3(fgm4) = 1;
    imagesc(I3)
    colormap('gray');
    %title('Modified regional maxima superimposed on original image (fgm4)')
    bw = imbinarize(Iobrcbr, 0.5);
    imagesc(bw);
    %title('Thresholded opening-closing by reconstruction (bw)')
    D = -bwdist(~bw);
    DL = watershed(D);
    bgm = DL == 0;
    %imagesc(bgm);
    %title('Watershed ridge lines (bgm)')
    gradmag2 = imimposemin(gradmag, bgm | fgm4);
    L = watershed(gradmag2);
    %L = imresize(L,1);
    %bgm = imresize(bgm,1);
    %fgm4 = imresize(fgm4,1);
    I4 = I;
    I4(imdilate(L == 0, ones(1, 1)) | bgm | fgm4) = 1;
    imagesc(I4)
    %title('Markers and object boundaries superimposed on original image (I4)')
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    imagesc(Lrgb)
    %L= imresize(L,0.5,'nearest');
    %title('Colored watershed label matrix (Lrgb)')
end

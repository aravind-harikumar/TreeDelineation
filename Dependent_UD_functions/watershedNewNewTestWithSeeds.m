function L = watershedNewNewTestWithSeeds(I,seedsIndx)
    %I =1-I;
   
    
    
    %I(1:50,1:2)=0;
    %I(1:50,49:50)=0;
    %I(49:50,1:50)=0;
    %I(1:2,1:50)=0;
    
    %load('imat2.mat')
    %I = imresize(I,2);
    %I = imgaussfilt(I,2);
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    imshow(gradmag,[])
    I = mat2gray(I);

    %title('Gradient magnitude (gradmag)')    
    se = strel('disk', 3);
    Io = imopen(I, se);
    imagesc(Io),
    Ie = imerode(I, se);
    imagesc(Ie),
    Iobr = imreconstruct(Ie, I);
    imagesc(Iobr)
    Ioc = imclose(Io, se);
    imagesc(Ioc)
    %title('Opening-closing (Ioc)')
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    imagesc(Iobrcbr)
    %title('Opening-closing by reconstruction (Iobrcbr)')
    %Iobrcbr = imdilate(Iobrcbr, se);
    %IobrcbrCopy(IobrcbrCopy<0.5) =0;
    %imagesc(Iobrcbr)
    
    seedsIndx = sub2ind(size(I), seedsIndx(:,2), seedsIndx(:,1));
    seedsIndx = uint16(seedsIndx);
    seeds = zeros(size(Iobrcbr));
    seeds(seedsIndx) = 1;
    seeds = imdilate(seeds, se);
    
    fgm = logical(seeds);%|imregionalmax(Iobrcbr,4);
    imagesc(fgm)
%     fgm = zeros(size(fgm));
%     fgm(69,18) = 1;
%     fgm(73,53) = 1;
%     fgm(63,82) = 1;
%     fgm = logical(imdilate(fgm,se));
%     imagesc(fgm)
    
    %title('Regional maxima of opening-closing by reconstruction (fgm)')
    I2 = I;
    I2(fgm) = 1;
    imagesc(I2)
    %title('Regional maxima superimposed on original image (I2)')
    se2 = strel(ones(1,1));
    fgm2 = imclose(fgm, se2);
    imagesc(fgm2)
    fgm3 = imerode(fgm2, se2);
    imagesc(fgm3)
    fgm4 = bwareaopen(fgm3, 20);
    imagesc(fgm4)
    I3 = I;
    I3(fgm4) = 1;
    imagesc(I3)
    %title('Modified regional maxima superimposed on original image (fgm4)')
    
    
    bw = imbinarize(Iobrcbr);
    imagesc(bw);
    %bw = double(bw);
    %bw(5:18,15:29) = 1;
    %bw =  imbinarize(bw);
    %imagesc(bw);
    
    %title('Thresholded opening-closing by reconstruction (bw)')
    D = -bwdist(~bw);    
    imagesc(D);

    %se = strel('disk',2,4);
    DL = watershed(D);
    bgm = DL == 0;
    imagesc(DL);
    
    gradsmag2 = imimposemin(gradmag, fgm4);
    imagesc(gradsmag2);
    %bgm = DL == 0;
    %imagesc(bgm);
    %title('Watershed ridge lines (bgm)')
    gradmag2 = imimposemin(gradmag, bgm | fgm4);
    gradmag2(isnan(gradmag2)) = 1;
    imagesc(gradmag2);
    L = watershed(gradmag2);
    %L = imresize(L,1);
    %bgm = imresize(bgm,1);
    %fgm4 = imresize(fgm4,1);
    imagesc(L);
    %I4 = I;
    %I4(imdilate(L == 0, ones(1, 1)) | bgm | fgm4) = 1;
    %imagesc(I4)
    %title('Markers and object boundaries superimposed on original image (I4)')
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    imagesc(Lrgb);
    
    %L = imresize(L,0.5,'nearest');
    
    %segResultd =  getOverlappedSegment(dsmImageSegmented, L);
    %imagesc(segResultd);
    %title('Colored watershed label matrix (Lrgb)')
end
function L = watershedNewNewTestIn(I,plotOn)
    %I =1-I;
   
%     I(1:50,1:5)=0;
%     I(1:50,45:50)=0;
%     I(45:50,1:50)=0;
%     I(1:5,1:50)=0;
    
    %load('imat2.mat')
    %I = imresize(I,2);
    %I = imgaussfilt(I,2);
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    
    if(plotOn)
        imshow(gradmag,[])
    end
    
    I = mat2gray(I);

    %title('Gradient magnitude (gradmag)')    
    se = strel('disk', 3);
    Io = imopen(I, se);
    if(plotOn)
        imagesc(Io);
    end
    Ie = imerode(I, se);
    if(plotOn)
        imagesc(Ie);
    end
    Iobr = imreconstruct(Ie, I);
    if(plotOn)   
        imagesc(Iobr);
    end
    Ioc = imclose(Io, se);
    if(plotOn)
        imagesc(Ioc);
    end
    %title('Opening-closing (Ioc)')
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    if(plotOn)
        imagesc(Iobrcbr);
    end
    %title('Opening-closing by reconstruction (Iobrcbr)')
    %Iobrcbr = imdilate(Iobrcbr, se);
    %IobrcbrCopy(IobrcbrCopy<0.5) =0;
    %imagesc(Iobrcbr)
    
    fgm = imregionalmax(Iobrcbr,4);
    if(plotOn)
        imagesc(fgm)
    end
%     fgm = zeros(size(fgm));
%     fgm(69,18) = 1;
%     fgm(73,53) = 1;
%     fgm(63,82) = 1;
%     fgm = logical(imdilate(fgm,se));
%     imagesc(fgm)
    
    %title('Regional maxima of opening-closing by reconstruction (fgm)')
    I2 = I;
    I2(fgm) = 1;
    if(plotOn)
        imagesc(I2)
    end
    %title('Regional maxima superimposed on original image (I2)')
    se2 = strel(ones(1,1));
    fgm2 = imclose(fgm, se2);
    if(plotOn)
        imagesc(fgm2);
    end
    fgm3 = imerode(fgm2, se2);
    if(plotOn)
        imagesc(fgm3);
    end
    fgm4 = bwareaopen(fgm3, 20);
    if(plotOn)
        imagesc(fgm4);
    end
    I3 = I;
    I3(fgm4) = 1;
    if(plotOn)
        imagesc(I3);
    end
    %title('Modified regional maxima superimposed on original image (fgm4)')
    bw = imbinarize(Iobrcbr);
    if(plotOn)
        imagesc(bw);
    end
    bw = double(bw);
    %bw(5:18,15:29) = 1;
    bw =  imbinarize(bw);
    if(plotOn)
        imagesc(bw);
    end
    
    %title('Thresholded opening-closing by reconstruction (bw)')
    D = bwdist(bw);
    if(plotOn)
        imagesc(D);
    end

    %se = strel('disk',2,4);
    DL = watershed(D);
    bgm = DL == 0;
    if(plotOn)
        imagesc(DL);
    end
    
    gradsmag2 = imimposemin(gradmag, fgm4);
    if(plotOn)
        imagesc(gradsmag2);
    end
    bgm = DL == 0;
    if(plotOn)
        imagesc(bgm);
    end
    %title('Watershed ridge lines (bgm)')
    gradmag2 = imimposemin(gradmag, bgm | fgm4);
    gradmag2(isnan(gradmag2)) = 1;
    if(plotOn)
        imagesc(gradmag2);
    end
    L = watershed(gradmag2);
    %L = imresize(L,1);
    %bgm = imresize(bgm,1);
    %fgm4 = imresize(fgm4,1);
    if(plotOn)
        imagesc(L);
    end
    I4 = I;
    I4(imdilate(L == 0, ones(1, 1)) | bgm | fgm4) = 1;
    if(plotOn)
        imagesc(I4);
    end
    %title('Markers and object boundaries superimposed on original image (I4)')
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    if(plotOn)
        imagesc(Lrgb);
    end
    
    %L = imresize(L,0.5,'nearest');
    
    %segResultd =  getOverlappedSegment(dsmImageSegmented, L);
    %imagesc(segResultd);
    %title('Colored watershed label matrix (Lrgb)')
end
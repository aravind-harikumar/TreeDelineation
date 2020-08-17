function L = watershedSeg(nDSM, seedss, showWShed, bachThresh)
        
        nDSM = mat2gray(nDSM);
        seeds = sub2ind(size(nDSM),seedss(:,2),seedss(:,1));
        [getX, getY] = ind2sub(size(nDSM),unique(seeds));        
        BW = zeros(size(nDSM));
        BW(unique(seeds)) = 255;
        se = strel('disk',2,4);
        BW = imdilate(BW,se);
        D = bwdist(BW,'euclidean'); 
        mask = imextendedmin(D,0.5);
        D2 = imimposemin(D,mask);
        
        L = watershed(D2);
        cellMsk = nDSM<(bachThresh);
        L(cellMsk)=0;
        if(showWShed)
            imagesc(L);
        end
        Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
        if(showWShed)
            imagesc(Lrgb), title('Opening-closing (Ioc)');
        end

        % % To get white boundary
        [SegList,Lbl] = bwboundaries(L,'noholes');
        Lbl(Lbl>0)=1; Lbl(Lbl==0)=2;
        colrmap = [105/255, 158/255, 112/255
                   227/255, 242/255, 229/255];
        if(showWShed)
            imagesc(label2rgb(Lbl, colrmap));
            hold on;
        end

        if(showWShed)
            plot(getY, getX, '.','MarkerSize',10,'Color','r');
            hold on;
        end

        for k = 1:length(SegList)
           boundary = SegList{k};
           if(showWShed)
                plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1.5)
           end
        end
        
        %title({'Tree tops from the CHM using LEVEL-SET algorithm.', ' Tree tops (in Red) Fixed Crown Extend (in Green)'});
end
function ITC_Proposed_Main()
    %% Clean Matlab environment %%      
    clear all; close all;
    
    %% Set dependencies folder paths %%    
    baseFolder = '.';
    addpath(genpath( strcat(baseFolder,'/Dependent_packages/') ));
    addpath(genpath( strcat(baseFolder,'/Dependent_UD_functions/') ));
    addpath(genpath( strcat(baseFolder,'/ITC_Segmentation/') ));
    
    %% Generate candidate segments (FYI: It overwrites the previous files if TRUE);
    run_CHM_3D_Segmentation = false;  cutTrees = false;
    inputLASPlotPath = 'Input_Files/Others/';
    outLASFilePath = fullfile(baseFolder,'/Output_files/las_files/ITCdata/'); % output_folder LS derived candidate segments
    % run segmentation algorithm %
    if(run_CHM_3D_Segmentation)
        disp('Obtainng 3D dominant tree data...');
        peaks_crh = runITC_WSSegmentation(inputLASPlotPath, outLASFilePath, cutTrees);
        figure; plot(peaks_crh(peaks_crh(:,4)==1,1), peaks_crh(peaks_crh(:,4)==1,2), '*')
    else
        disp(strcat('Using existing segmentation reults...'));
    end
    
    %% Perform Sub-dominant tree detection for each candidate segment detected in a plot %%
    files = dir(fullfile(outLASFilePath, strcat('')));  
    treeCount = 1;
    Final_op_Arr=[];
    for file = files'
        segmentName = file.name;
        if( strcmp(segmentName, '1_sub_tree_1_dom_tree.las') )
            lasFileName =  fullfile(outLASFilePath, segmentName);
            [t_topArr, X_ExtendArr,Y_ExtendArr, stData] =  ITC_Proposed_Algo(lasFileName, outLASFilePath, treeCount);
            %getTreeTops = []; getXExtTree = []; getYExtTree = []; XY_Dom_Tree = []; opArr ={};
            if(size(t_topArr,1)>0)
                opArr{1} = segmentName
                Final_op_Arr = [Final_op_Arr;...
                                table(repmat(opArr,size(t_topArr,1),1), ...
                                t_topArr,...
                                X_ExtendArr,...
                                Y_ExtendArr,...
                                repmat(stData.origDataMaxXY,size(t_topArr,1),1))];
            end
            treeCount=treeCount+1;
            pause(1);
       end
    end    
    % Delete the file and write again
    if(exist('Sub_Tree_Data.xls','file')==2)
        delete 'Sub_Tree_Data.xls';
    end
%     writetable(Final_op_Arr,'Sub_Tree_Data.xls');
end

function segPointsArr = getliDARSegmentsLM(stData, segmentArr, Lbl)
    
    [maxX, maxY] = deal(max(stData.lidarDataArray(:,1)),max(stData.lidarDataArray(:,2)));
    [minX, minY] = deal(min(stData.lidarDataArray(:,1)),min(stData.lidarDataArray(:,2)));
    xStep = (maxX - minX)/ size(Lbl,1);
    %yStep = (maxY - minY)/ size(Lbl,2);
    
    lData = stData.lidarDataArray(:,1:3);
    maxDDX = max(lData(:,1));
    maxDDY = max(lData(:,2));
    
    %axis([-250 250 -250 250 -50 50])
    lData(:,1) = lData(:,1)-minX;
    lData(:,2) = maxDDY - (lData(:,2));
    
%     J = scatter3(lData(:,1)/xStep,lData(:,2)/xStep,lData(:,3),150,lData(:,3),'.');      
% 
%     T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
%     x = [0 85 110 255];
%     map = interp1(x/255,T,linspace(0,1,100));
%     colormap(map)

    
    segPointsArr =[];
    for iSegCount = 1:1:size(segmentArr,1)
        
        subplot(2,2,[1]);
        boundary = segmentArr{iSegCount};
        [inV,onV] = inpolygon(lData(:,1)/xStep,lData(:,2)/xStep,boundary(:,2), boundary(:,1)); 
        
        AllIndices = or(inV,onV);
        %plot(lData(AllIndices,1)/xStep,lData(AllIndices,2)/xStep,'*');

        plot(boundary(:,2), boundary(:,1),'.-.','Color',[0 1 0], 'LineWidth', 2.5);    

        segPointsArr{iSegCount} = stData.lidarDataArray(AllIndices,6);
        
        f=subplot(2,2,[2]);
        h(1) = scatter3(lData(AllIndices,1),lData(AllIndices,2), lData(AllIndices,3),150, lData(AllIndices,3),'.');
        axis equal; grid on;
        T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
        x = [0 85 110 255];
        map = interp1(x/255,T,linspace(0,1,100));
        colormap(f,map);
        
        camproj perspective; rotate3d on; axis vis3d; view(90 ,0);
        %axis([-100 100 -100 100 0 50]);
        if(mod(iSegCount,50))
            fprintf('.');
        else
            fprintf('\n');
        end         
        pause(0.01);
    end 
    fprintf('\n');
end


function segPointsArr = getliDARSegmentsLS(stData, segmentArr, Lbl)
    
    [maxX, maxY] = deal(max(stData.lidarDataArray(:,1)),max(stData.lidarDataArray(:,2)));
    [minX, minY] = deal(min(stData.lidarDataArray(:,1)),min(stData.lidarDataArray(:,2)));
    xStep = (maxX - minX)/ size(Lbl,1);
    %yStep = (maxY - minY)/ size(Lbl,2);
    
    lData = stData.lidarDataArray(:,1:3);
    maxDDX = max(lData(:,1));
    maxDDY = max(lData(:,2));
    
    %axis([-250 250 -250 250 -50 50])
    lData(:,1) = lData(:,1)-minX;
    lData(:,2) = maxDDY - (lData(:,2));
    
%     J = scatter3(lData(:,1)/xStep,lData(:,2)/xStep,lData(:,3),150,lData(:,3),'.');      
% 
%     T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
%     x = [0 85 110 255];
%     map = interp1(x/255,T,linspace(0,1,100));
%     colormap(map)

    
    segPointsArr =[];
    for iSegCount = 1:1:size(segmentArr,1)
        
        subplot(2,2,[3]);
        
        boundary = segmentArr{iSegCount}; magFactor = 0;
        centerOrig = [mean(boundary(:,1)) mean(boundary(:,2))];        
        boundaryBig = [boundary(:,1)+boundary(:,1)*magFactor boundary(:,2)+boundary(:,2)*magFactor];        
        centerBB = [mean(boundaryBig(:,1)) mean(boundaryBig(:,2))];
        diffX = centerBB(1) -  centerOrig(1);
        diffY = centerBB(2) -  centerOrig(2);        
        boundaryBigShifted = [boundary(:,1)+boundary(:,1)*magFactor-diffX boundary(:,2)+boundary(:,2)*magFactor-diffY]; 
        boundary = boundaryBigShifted;
        %figure;
        %plot(boundary(:,2),boundary(:,1),'*')
        %hold on;
        boundary(:,2) = smooth(boundary(:,2),15);
        boundary(:,1) = smooth(boundary(:,1),15);
        %plot(boundary(:,2),boundary(:,1),'*')
        
        
        [inV,onV] = inpolygon(lData(:,1)/xStep,lData(:,2)/xStep,boundary(:,2), boundary(:,1)); 
        
        AllIndices = or(inV,onV);
        %plot(lData(AllIndices,1)/xStep,lData(AllIndices,2)/xStep,'*');

        plot(boundary(:,2), boundary(:,1),'.-.','Color',[0 1 0], 'LineWidth', 2.5);    

        segPointsArr{iSegCount} = stData.lidarDataArray(AllIndices,6);
        
        f=subplot(2,2,[4]);
        h(1) = scatter3(lData(AllIndices,1),lData(AllIndices,2), lData(AllIndices,3),150, lData(AllIndices,3),'.');
        axis equal; grid on;
        T = [10, 6, 255;  0 255 0;  255, 255, 0; 255, 0, 0]./255;
        x = [0 85 110 255];
        map = interp1(x/255,T,linspace(0,1,100));
        colormap(f,map)
        
        camproj perspective; rotate3d on; axis vis3d; view(90 ,0);
        %axis([-100 100 -100 100 0 50]);
        if(mod(iSegCount,50))
            fprintf('.');
        else
            fprintf('\n');
        end         
        pause(0.01);
    end 
    fprintf('\n');
end


function SegList = localMaxbasedSegmentation(plotID, nDSM, segmentationAlgo, showDSM)
    
    SegList = [];
    %nDSMasd = imread(fullfile(inSeedFilepath, '',strcat(plotID,'_seed.tif')));    
    %maxindxes =find(nDSMasd>0);
    %nDSM = imgaussfilt(nDSM,1);
    
    [cent, varargout] = FastPeakFind(mat2gray(nDSM),10);
    
    getX = cent(1:2:end);
    getY = cent(2:2:end);
    
    % print the tree tops
    imshow(nDSM,[0 255]);
    hold on;
    %[getX, getY] = ind2sub(size(nDSM),unique(maxindxes));
    plot(getX, getY, '.','MarkerSize',10,'Color','r');
    % Plot circles
    maxindxes = sub2ind(size(nDSM),getY,getX);
    % Perform Segmentation    
    if(strcmp(segmentationAlgo,'circular'))
        SegList =[]; Lbl= zeros(size(nDSM));
        for j =1:1:size(getY,1)
            [ypFinal, xpFinal] = circle(getX(j),getY(j),25);
            jointArr  = [ypFinal' xpFinal'];
            SegList{j} =  jointArr;
            hold on;
        end
        SegList = SegList';
        axis square; 
        title({'Tree tops from the CHM using LOCAL_MAX algorithm.', ' Tree tops (in Red) Fixed Crown Extend (in Green)'});
        hold on;
    else                
        disp(strcat('Performing watershed on nDSM...'));        
        L=watershedNewNewTest(imgaussfilt(nDSM,0.1),maxindxes, true);
        hold on;
        colormap(colorcube);
        plot(getX, getY, '.','MarkerSize',20,'Color','r');
        axis square; 
        set(gcf, 'Position', get(0,'Screensize')); 
        saveas(gcf,  strcat(plotID,'_segmented_map_LM.tif'))
% 
%                 [SegList,Lbl] = bwboundaries(L,'noholes');
%                 Lbl(Lbl>0)=1; Lbl(Lbl==0)=2;
%                 hold on;
%                 for k = 1:length(SegList)
%                     boundary = SegList{k};
%                     magFactor = 0.0;
%                     centerOrig = [mean(boundary(:,1)) mean(boundary(:,2))];        
%                     boundaryBig = [boundary(:,1)+boundary(:,1)*magFactor boundary(:,2)+boundary(:,2)*magFactor];        
%                     centerBB = [mean(boundaryBig(:,1)) mean(boundaryBig(:,2))];
%                     diffX = centerBB(1) -  centerOrig(1);
%                     diffY = centerBB(2) -  centerOrig(2);        
%                     boundaryBigShifted = [boundary(:,1)+boundary(:,1)*magFactor-diffX boundary(:,2)+boundary(:,2)*magFactor-diffY]; 
%                     boundary = boundaryBigShifted;
%                     %plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
%                 end        
%                 axis square;
    end
       
end


function watershedSeg(nDSM, seeds)

        [getX, getY] = ind2sub(size(nDSM),unique(seeds));

        %nDSMasd = imread(fullfile(inSeedFilepath, 'LevelSetSeeds',strcat(plotID,'_seed.tif')));    
        %maxindxes =find(nDSMasd>0);
        subplot(2,2,[1]);
        BW = zeros(size(nDSM));
        BW(unique(seeds)) = 255;
        se = strel('disk',2,4);
        BW = imdilate(BW,se);
        imagesc(BW);

        D = bwdist(BW,'euclidean'); 
        mask = imextendedmin(D,0.5);
        imshow(D,[min(D(:)) max(D(:))]);
        D2 = imimposemin(D,mask);
        imagesc(D2);

        L = watershed(D2);
        cellMsk = nDSM<(0*255);
        L(cellMsk)=0;
        imagesc(L), 
        Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
        ax2=subplot(2,2,[1 3]);
        imshow(Lrgb,[0 255]), title('Opening-closing (Ioc)');

        % % To get white boundary
        [SegList,Lbl] = bwboundaries(L,'noholes');
        Lbl(Lbl>0)=1; Lbl(Lbl==0)=2;
        ax2=subplot(2,2,[1]);
        colrmap = [105/255, 158/255, 112/255
                   227/255, 242/255, 229/255];
        imshow(label2rgb(Lbl, colrmap));
        hold on;

        plot(getY, getX, '.','MarkerSize',10,'Color','r');
        hold on;

        for k = 1:length(SegList)
           boundary = SegList{k};
           plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1.5)
        end
        
        axis square; 
        title({'Tree tops from the CHM using LEVEL-SET algorithm.', ' Tree tops (in Red) Fixed Crown Extend (in Green)'});
end

%%
function SegList = levelSetbasedSegmentation(inSeedFilepath, plotID, nDSM, segmentationAlgo, showDSM)
    
    SegList = {};
    nDSMasd = imread(fullfile(inSeedFilepath, '',strcat(plotID,'_seed.tif')));
    tempnDSMasd =  fliplr(imrotate(nDSMasd,-90));

    CC = bwconncomp(tempnDSMasd);
    S = regionprops(CC,'Centroid');
    centroidpairsArray = struct2array(S);        
    getX = centroidpairsArray(1:2:end)';
    getX = round(getX,0);
    getY = centroidpairsArray(2:2:end)';
    getY = round(getY,0);
    maxindxes = sub2ind(size(tempnDSMasd), getY, getX);

    hold on;
    nDSMasd(1,:) = 0; nDSMasd(:,1) = 0; nDSMasd(end,:) = 0; nDSMasd(:,end) = 0;

    % print the tree tops
    imshow(nDSM,[0 255]);
    hold on;
    %[getX, getY] = ind2sub(size(nDSM),unique(maxindxes));
    plot(getY, getX, '.','MarkerSize',15,'Color','r');
    % Plot circles
    hold on;
    
    % Perform Segmentation    
    if(strcmp(segmentationAlgo,'circular'))
        SegList =[]; Lbl= zeros(size(nDSM));
        for j =1:1:size(getY,1)
            [ypFinal, xpFinal] = circle(getX(j),getY(j),25);
            jointArr  = [ypFinal' xpFinal'];
            SegList{j} =  jointArr;
            hold on;
        end
        SegList = SegList';
        axis square; 
        title({'Tree tops from the CHM using LEVEL-SET algorithm.', ' Tree tops (in Red) Fixed Crown Extend (in Green)'});
    else
        
        disp(strcat('Performing watershed on nDSM...'));   
        
        CC = bwconncomp(nDSMasd);
        S = regionprops(CC,'Centroid');
        centroidpairsArray = struct2array(S);        
        getXX = centroidpairsArray(1:2:end)';
        getXX = round(getXX,0);
        getYY = centroidpairsArray(2:2:end)';
        getYY = round(getYY,0);
        maxindxes = sub2ind(size(nDSMasd), getYY, getXX);
   
        imagesc(nDSM);
        c1=[0 0 0]; c2=[50/256 77/256 0/256]; c3=[122/256 190/256 0/256]; c4=[225/256 231/256 221/256]; plotInternalOn = false;
        cmap = getColorMap(c1, c2, c3, c4, size(colormap) , plotInternalOn);
        colormap(cmap);
        
        plotOn = false;
        L=watershedNewNewTest(nDSM,maxindxes,plotOn);
        uniqueLInds = unique(L);        

        myArr = []; mxArr =[];  EccentricityArr=[];
        for k = 1:length(uniqueLInds)
             pixIndxSegment = find(L==k);
             if(length(pixIndxSegment)>0)
                 [ccxIndx, ccyIndex] =  ind2sub(CC.ImageSize, pixIndxSegment);
                 mx = mean(ccxIndx);        
                 my = mean(ccyIndex);             
                 tempBinImage = zeros(CC.ImageSize); 
                 tempBinImage(pixIndxSegment) = 1;             
                 S = regionprops(tempBinImage,'Eccentricity');
                 SBB = regionprops(tempBinImage,'BoundingBox');
                 SBBArea = SBB.BoundingBox(3)*SBB.BoundingBox(4);
                 
                 se = strel('disk', 1);
                 if(sum(tempBinImage(:)==1)>0)
                    tempBinImage = imdilate(tempBinImage, se);
                 end
                 %Io = imerode(tempBinImage, se);
                
                 [boundaryCell,~] = bwboundaries(tempBinImage,'noholes');
                 SegList{k} = boundaryCell;
                 
                 boundaryArray = boundaryCell{1};
                 plot(boundaryArray(:,2), boundaryArray(:,1), '-', 'LineWidth', 2, 'Color',[1 1 1])
                  
                 % to get rid of background noisy segments.
                 if(SBBArea > 0.08*(CC.ImageSize(1)*CC.ImageSize(2)))
                    L(pixIndxSegment) = 0;    
                 end

                 if(S.Eccentricity<0) % 0 is close to circle, 1 = close to line; if(S.Eccentricity < 0.6)
                    L(pixIndxSegment) = 0;   
                 else
                     
                     myArr = [myArr my];
                     mxArr = [mxArr mx];
                     EccentricityArr = [EccentricityArr S.Eccentricity];
                 end             
             end
        end
        
        %imagesc(L);
        %cmap = colorcube;
        %cmap(1,:) = [1 1 1];
        %cmap(end,:) = [0 0 0];
        %colormap(cmap);
        
        for j =1:1:length(myArr)
           % text(myArr(j), mxArr(j), num2str(EccentricityArr(j)),'HorizontalAlignment','right', 'FontSize',6,'Color','k'); 
        end        
        axis square;

        hold on;
        axis square;
        plot(getY, getX, '.','MarkerSize',20,'Color','r'); 
        axis square;
        set(gcf, 'Position', get(0,'Screensize')); 
        %saveas(gcf,  strcat(plotID,'_segmented_map_LS.tif'))
        
% %       To get white boundary
%         [SegList,Lbl] = bwboundaries(L,'noholes');
%         Lbl(Lbl>0)=1; Lbl(Lbl==0)=2;
%         hold on;
% 
%         for k = 1:length(SegList)
%             boundary = SegList{k};
%             magFactor = 0.0;
%             centerOrig = [mean(boundary(:,1)) mean(boundary(:,2))];        
%             boundaryBig = [boundary(:,1)+boundary(:,1)*magFactor boundary(:,2)+boundary(:,2)*magFactor];        
%             centerBB = [mean(boundaryBig(:,1)) mean(boundaryBig(:,2))];
%             diffX = centerBB(1) -  centerOrig(1);
%             diffY = centerBB(2) -  centerOrig(2);        
%             boundaryBigShifted = [boundary(:,1)+boundary(:,1)*magFactor-diffX boundary(:,2)+boundary(:,2)*magFactor-diffY]; 
%             boundary = boundaryBigShifted;
%             plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
%         end
%         %plot(getY, getX, '.','MarkerSize',5,'Color','r');
%         axis square; 
%         alpha(0.9);
%         set(gcf, 'Position', get(0, 'Screensize'));
%         title({'Tree crown boundaries'});
    end
end
%%


function L = watershedNewNewTest(I,seedsIndx,plotOn)

    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    if(plotOn)
        imshow(gradmag,[])
    end
    I = mat2gray(I);

    if(plotOn)
        hold off;
    end 
    %title('Gradient magnitude (gradmag)')    
    se = strel('disk', 1);
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
        imagesc(Iobr)
    end
    Ioc = imclose(Io, se);
    if(plotOn)
        imagesc(Ioc)
    end
    %title('Opening-closing (Ioc)')
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    if(plotOn)
        imagesc(Iobrcbr)
    end

    fgm = imregionalmax(Iobrcbr,4);
    if(plotOn)
        imagesc(fgm)
    end
   
    seeds = zeros(size(I));
    seeds(seedsIndx) = 1;
    se = strel('disk', 1);
    seeds = imdilate(seeds, se);
    fgm = logical(seeds);
    %title('Regional maxima of opening-closing by reconstruction (fgm)')
    I2 = I;
    I2(fgm) = 1;
    if(plotOn)
        imagesc(I2)
    end
    %title('Regional maxima superimposed on original image (I2)')
    se2 = strel('disk', 1);
    fgm2 = imclose(fgm, se2);
    if(plotOn)
        imagesc(fgm2)
    end
    fgm3 = imerode(fgm2, se2);
    if(plotOn)
        imagesc(fgm3)
    end
    fgm4 = bwareaopen(fgm3, 1);
    if(plotOn)
        imagesc(fgm4)
    end
    I3 = I;
    I3(fgm4) = 1;
    if(plotOn)
        imagesc(I3)
    end
    %title('Modified regional maxima superimposed on original image (fgm4)')
    
    
    bw = imbinarize(Iobrcbr);
    if(plotOn)
        imagesc(bw);
    end
    bw = double(bw);
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
    L(L==0) = 99999; % A BIG NUMBER
    cellMsk = I<(0.1);
    L(cellMsk)=0;
    if(plotOn)
        imagesc(L)
    end
    I4 = I;
    I4(imdilate(L == 0, ones(1, 1)) | bgm | fgm4) = 1;
    %imagesc(I4)
    %title('Markers and object boundaries superimposed on original image (I4)')
    %Lrgb = label2rgb(L, 'jet', 'k', 'shuffle');
    %imagesc(Lrgb);
    %hold on;
end


function cropDominantTreeData(segPointsArr, fullFileName, outFilepath)    
    % cut out individual trees (i.e., from input point cloud)
    for jCount = 1:1:length(segPointsArr)
        
        plotLevelLiDARdata = LoadData(fullFileName);
        stData = dataPreProcessing(plotLevelLiDARdata);
        
        % get LiDAR points belonging to a segment
        lidarDataArraytemp = ismember(stData.lidarDataArray(:,6),segPointsArr{jCount});
        if(sum(lidarDataArraytemp)>0) % if segments has at least 1 LiDAR point data     
            treeData = plotLevelLiDARdata; % for initializing with same structure
            %plotLevelLiDARdata = plotLevelLiDARdata; % for getting same class object            
            treeData.x = plotLevelLiDARdata.x(lidarDataArraytemp);
            treeData.y = plotLevelLiDARdata.y(lidarDataArraytemp);
            treeData.z = plotLevelLiDARdata.z(lidarDataArraytemp);
            intensity = plotLevelLiDARdata.get_intensity();
            treeData.intensity = intensity(lidarDataArraytemp);
            treeData.bits =[];
            classification = plotLevelLiDARdata.get_classification();
            treeData.classification = classification(lidarDataArraytemp);
            user_data = plotLevelLiDARdata.get_user_data();
            treeData.user_data = user_data(lidarDataArraytemp);
            scan_angle = plotLevelLiDARdata.get_scan_angle();
            treeData.scan_angle = scan_angle(lidarDataArraytemp);
            point_source_id = plotLevelLiDARdata.get_point_source_id();
            treeData.point_source_id = point_source_id(lidarDataArraytemp);
            gps_time = plotLevelLiDARdata.get_gps_time();
            try
                treeData.gps_time = gps_time(lidarDataArraytemp);
            catch
                ss = 0;
            end
            treeData.read_color();
            treeData.red=[];
            treeData.green=[];
            treeData.blue=[];
            treeData.nir=[];
            treeData.extradata = [];
            treeData.read_point_wave_info();
            treeData.Xt=[];
            treeData.Yt=[];
            treeData.Zt=[];
            treeData.wavedescriptors=[];    
            treeData.selection = true( size(treeData.x,1),1);
            treeData.wave_return_point=[];
            treeData.extradata=[];

            % write .las file.
            obj = write_las(treeData, strcat(outFilepath ,'Segment', num2str(jCount),'.las'));
            fclose('all');
           
            if(mod(jCount,50))
                fprintf('.');
            else
                fprintf('\n');
            end
            
        end

    end
    fprintf('\n');
end

%%
function returnData =  LoadData(fullFileName)
    returnData = lasdata(fullFileName);
end
%%
function [retMaxXYZ, indx] = findMaxHeightXY(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
end
%%
function retTable = write2table(lasFile)
    retTable = zeros(size(lasFile.x,1),5);
    retTable(:,1) = lasFile.x;
    retTable(:,2) = lasFile.y;
    retTable(:,3) = lasFile.z;
    retTable(:,4) = get_classification(lasFile);    
    retTable(:,5) = lasFile.get_return_number;
end
%%
function crownHt = getMinCrownHeight(lidarDataArray)
   maxTreeHeight = max(lidarDataArray(:,3));
   crownHt = min(lidarDataArray(:,3));
   for i = maxTreeHeight:-1:1
       templidararray = lidarDataArray(find(and(lidarDataArray(:,3)<i,(lidarDataArray(:,3)>i-1))),1:3);
        if(size(templidararray,1)<10)
            crownHt = i;
            break;
        end
   end
end
%%
function dataAtt = dataPreProcessing(singleTreeLiDARdata)
% Write normalized data to a table for performance improvement 
        dataAtt.origData = write2table(singleTreeLiDARdata);

        lidarDataArr = normalizeLiDARData(write2table(singleTreeLiDARdata));
        %lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:);
               
        dataAtt.treeWidth =  max(lidarDataArr(:,1)); treeHeight = max(lidarDataArr(:,2)); % to calculate max wisth/breadth to set plot width
        dataAtt.maxTreeWidth = max(dataAtt.treeWidth, treeHeight); % Keep it as 5 if issues arise
        dataAtt.mintreeHeight = min(lidarDataArr(:,3));
        dataAtt.maxtreeHeight = max(lidarDataArr(:,3));        
        dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;        
        dataAtt.htDeduction = 0; %dataAtt.treeHeight*0.05; % to get rid of ground noise points
        
        % Crown height
        lidarDataArr = lidarDataArr(find(and(lidarDataArr(:,3) >= min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) <= max(lidarDataArr(:,3)))),:);
        dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
        dataAtt.maxCrownHeight = max(lidarDataArr(:,3));        
        dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;       
              
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 > 0.1),:);                
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);
        
        % Identify the center point of the tree (in top view).
        [dataAtt.retMaxXYZ, dataAtt.retMaxXYZIndex] = findMaxHeightXY(lidarDataArr);
        
        %lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); 
        %lidarDataArr = lidarDataArr(find(lidarDataDensityArr(:,1) >= 0.0000),:);
        
        %lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); %get point in very close neightbourhood                
        index = 1:1:size(lidarDataArr,1); 
         
        % % Density priuned LiDAR data
        %selIndx = find(lidarDataDensityArr(:,1) >= 0.0000);        
        dataAtt.lidarDataArray =  lidarDataArr; %lidarDataArr(selIndx,:);                
        dataAtt.index = index; %index(selIndx);
        dataAtt.lidarDataDensityArr = zeros(size(lidarDataArr,1),1); %lidarDataDensityArr(selIndx,:);        
        dataAtt.lidarDataArray = [dataAtt.lidarDataArray dataAtt.index' dataAtt.lidarDataDensityArr];
        
        dataAtt.lidarDataArray(:,1) = dataAtt.lidarDataArray(:,1) - dataAtt.retMaxXYZ(1);
        dataAtt.lidarDataArray(:,2) = dataAtt.lidarDataArray(:,2) - dataAtt.retMaxXYZ(2);
        
        
        [dataAtt.retMaxXYZNew, dataAtt.retMaxPointIndexNew] = findMaxHeightXYNew(dataAtt.lidarDataArray);

end
%%
function lidarDataArray = normalizeLiDARData(lidarDataArray)
    midxPoint = (max(lidarDataArray(:,1))- min(lidarDataArray(:,1)))/2;
    midyPoint = (max(lidarDataArray(:,2))- min(lidarDataArray(:,2)))/2;
    
    lidarDataArray(:,1) = lidarDataArray(:,1)- min(lidarDataArray(:,1));
    lidarDataArray(:,2) = lidarDataArray(:,2)- min(lidarDataArray(:,2));
    lidarDataArray(:,3) = lidarDataArray(:,3)- min(lidarDataArray(:,3));

    lidarDataArray(:,1) = lidarDataArray(:,1)- midxPoint;
    lidarDataArray(:,2) = lidarDataArray(:,2)- midyPoint;
end
%%
function [retMaxXYZ, maxpointIndex] = findMaxHeightXYNew(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
     maxpointIndex = lidarDataArray(indx,6);
end
%%
function [xpFinal, ypFinal] = circle(x,y,r)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.2:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    xpFinal = x+xp; ypFinal = y+yp;
    plot(ypFinal,xpFinal,'.-','MarkerSize',2,'Color',[111/255, 255/255, 111/255]);
end
% %%

% function ITC_Proposed_Main()
%     %% Clean Matlab environment %%      
%     clear all; close all;
%     
%     %% Set dependencies folder paths %%    
%     baseFolder = '.';
%     addpath(genpath( strcat(baseFolder,'\Dependent_packages\') ));
%     addpath(genpath( strcat(baseFolder,'\Dependent_UD_functions\') ));  
%     
%     %% Set folder structure %%   
%     inLASFilepath = fullfile(baseFolder,'\Input_Files\plot_las_files\NEW\'); % input folder_plot_level
%     inSeedFilepath = fullfile(baseFolder,'\Input_Files\plot_seeds\LevelSetSeeds\NEW\'); % seed_folder_plot_level
%     %outLASFilepathLM = fullfile(baseFolder,'\Output_files\las_files\LM\'); % output_folder LM derived candidate segments
%     outLASFilepathLS = fullfile(baseFolder,'\Output_files\las_files\for_paper_real_data_results\'); % output_folder LS derived candidate segments
%     plotID = 'Ads6_dcm';
%     nDSMImg = imread(fullfile(baseFolder,'\Input_Files\plot_nDSMs\NEW\CHMS_BW_ORIGINAL\',strcat(plotID,'.tif'))); % nDSM file    
%     nDSMImg = nDSMImg(:,:,1);
%     nDSMImg = imgaussfilt(nDSMImg,0.5);
%     %nDSMImg(nDSMImg<10) = 0;
%     % Generate candidate segments (FYI: It overwrites the previous files if TRUE);
%     run_CHM_3D_Segmentation = false;    
%     
%     % Load LiDAR data%%
%     disp(strcat('Loading LiDAR data...'));
%     fullFileName = strcat(inLASFilepath,plotID,'.las');
%     standLiDARdata = LoadData(fullFileName);
%     
%     % LiDAR data Pre-processing and basic tree information retrieval %%
%     if(run_CHM_3D_Segmentation)
%         disp(strcat('Preprocessing data...'));
%         stData = dataPreProcessing(standLiDARdata);
%         % Obtain candidate segments %%
%         disp(strcat('Segmenting dominant trees...'));
%         % perform segmentation on nDSM
%         disp(strcat('Executing Step 1: obtaining tree tops'));
%         
%         %Seed detection algorithm from CHM can be 'Local_Maxima' or 'Level_Set'
%         showDSM = true;
%         %segmentArrLM = localMaxbasedSegmentation(plotID, nDSMImg, 'noncircular', showDSM);
%         segmentArrLS = levelSetbasedSegmentation(inSeedFilepath, plotID, nDSMImg, 'noncircular', showDSM); % watershed or circular (by max radius)
% 
%         % Segments out corresponding 3D tree point clouds
%         disp(strcat('Executing Step 2: Obtainng 3D dominant tree data...'));
%         %segPointsArrLM = getliDARSegmentsLM(stData, segmentArrLM, nDSMImg); %3D segmentation        
%         segPointsArrLS = getliDARSegmentsLS(stData, segmentArrLS, nDSMImg); %3D segmentation  
% 
%         disp(strcat('Executing Step 3: Writing 3D dominant tree data as .las files...'));
%         %cropDominantTreeData(segPointsArrLM, fullFileName, outLASFilepathLM); %3D point cloud writing
%         cropDominantTreeData(segPointsArrLS, fullFileName, outLASFilepathLS); %3D point cloud writing
%     else
%         disp(strcat('Using existing segmentation reults...'));
%     end
%     % 1_sub_tree_1_dom_tree_C9 1_sub_tree_1_dom_tree_C7 000370_elli_elli_Tile_8_2trees_20 - 1 tree
%     
%     %% Perform Sub-dominant tree detection for each candidate segment detected in a plot %%
%     files = dir(fullfile(outLASFilepathLS, strcat('*.las')));  
%     treeCount = 1;
%     Final_op_Arr=[];
%     for file = files'
%         segmentName = file.name;
%         if( strcmp(segmentName, '11.las') )
%             lasFileName =  fullfile(outLASFilepathLS, segmentName);            
%             [t_topArr, X_ExtendArr,Y_ExtendArr, stData] =  ITC_Proposed_Algo(lasFileName, outLASFilepathLS,treeCount);
%             
%             %getTreeTops = []; getXExtTree = []; getYExtTree = []; XY_Dom_Tree = []; opArr ={}; 
%             if(size(t_topArr,1)>0)            
%                 opArr{1} = segmentName  
%                 Final_op_Arr = [Final_op_Arr;...
%                                 table(repmat(opArr,size(t_topArr,1),1), ...
%                                 t_topArr,...
%                                 X_ExtendArr,...
%                                 Y_ExtendArr,...
%                                 repmat(stData.origDataMaxXY,size(t_topArr,1),1))];
%             end
%             treeCount=treeCount+1; pause(1);
%        end
%     end
%     % Delete the file and write again
%     if(exist('Sub_Tree_Data.xls','file')==2)
%         delete 'Sub_Tree_Data.xls';
%     end        
%     writetable(Final_op_Arr,'Sub_Tree_Data.xls');
%     
% end

function [f3Dvoxel, voxelFeatureWiseIndices] = getvoxelGLCMFeatures(f3Dvoxel, VOX_DIVS, plot3D, crownSegmentsAll, domTreeThrshLevel, colorArr)
           
        subTreeVoxels = zeros(VOX_DIVS.xDiv, VOX_DIVS.yDiv, VOX_DIVS.zDiv);        
        VoxelfeatureArr =[]; numFeatures = 12;
        getSubMatrices = get3DIndices(f3Dvoxel.voxNeighArrDensity,[VOX_DIVS.xDiv VOX_DIVS.yDiv VOX_DIVS.zDiv]);

        fprintf('\n');
        disp(strcat({'   Building voxels space...'}));

        % Parameters of GLCM
        totalGLCMFeatureCount = 12;
        GLCMDistances = [1]; 
        GLCMDirections = [ 0  1   0;   -1  1  0;   -1  0  0;    -1 -1  0; ...
                           0  1  -1;    0  0 -1;    0 -1 -1;    -1  0 -1; ...
                           1  0  -1;   -1  1 -1;    1 -1 -1;    -1 -1 -1;  1  1 -1];
        indxArr = getIc3d(size(GLCMDistances,1), size(GLCMDirections,1), totalGLCMFeatureCount); %  size(GLCMDirections,1)

        %fNames = fieldnames(f3Dvoxel);
        %fet2skip = size(fNames,1); % selected featire = fetcnt+1;%         
        %load('matlab_f_1.mat');        
        fNames = {'voxNeighArrDensity'};
        % GLCM feature extraction
%         fNames = {'voxNeighArrDensity','energ','Corr','contr','homom'};
%         maxCnt = size(getSubMatrices,2);
%         for j =1:1:size(getSubMatrices,2)
%             [tt,gg,dd]=ind2sub([VOX_DIVS.xDiv VOX_DIVS.yDiv VOX_DIVS.zDiv],j);
%             [textureArr, hhh] = cooc3d(getSubMatrices{j}, 'distance', GLCMDistances, 'direction', GLCMDirections);
%             %fNames = {'Corr','contr','corrm','corrp','cprom','cshad','dissi','energ','entro','homom','homop','maxpr'};
%             for fCount = 1:size(fNames,2)
%                 f3Dvoxel.(fNames{fCount})(tt,gg,dd) =  mean(textureArr(indxArr(fCount,:)));
%             end
%             if(~mod(j, ceil(maxCnt/25)))
%                 perccount(j,maxCnt);
%             end
%         end
        
       % load('matlab.mat');       
        texFeatureArr = zeros(VOX_DIVS.xDiv, VOX_DIVS.yDiv, VOX_DIVS.zDiv);
%         figure; %kmeans based analaysis
%         for i = 1:1:VOX_DIVS.zDiv           
%             feature0 = f3Dvoxel.('voxNeighArrDensity')(:,:,i);
%             feature1 = f3Dvoxel.('energ')(:,:,i);
%             feature2 = f3Dvoxel.('Corr')(:,:,i);
%             feature3 = f3Dvoxel.('contr')(:,:,i);
%             feature4 = f3Dvoxel.('homom')(:,:,i);            
%             KMInput = [feature0(:) feature1(:) feature2(:) feature3(:) feature4(:)];
%             KMInput(isnan(KMInput)) = 0;
%             % repeat the clustering 3 times to avoid local minima
%             [cluster_idx, cluster_center] = kmeans(KMInput,4,'distance','sqEuclidean','Replicates',3);   
%             cluster_idx = reshape(cluster_idx,size(feature1,1),size(feature1,2));
%             csmcdm = flipud(mat2gray(cluster_idx));
%             subplot(7,7,i); imagesc(csmcdm);
%         end
        
        % Corr contr corrm corrp cprom cshad dissi energ entro homom homop maxpr
        RowArr =[];
        for icnt = 1:size(crownSegmentsAll,2)
           [aRow,bRow] = ind2sub( size(crownSegmentsAll{icnt}),find(crownSegmentsAll{icnt}>0) );
           if(~ isempty(aRow))
                RowArr = [RowArr;round(mean(bRow),0) round(mean(aRow),0) 255];
           end
        end
        
       % RowArr(:,1) = round((RowArr(:,1)/size(crownSegmentsAll{1},2))*VOX_DIVS.xDiv,0);
       % RowArr(:,2) = round((RowArr(:,2)/size(crownSegmentsAll{1},1))*VOX_DIVS.yDiv,0);
        
        RowArr = [RowArr; 5 5 255];
        %[VOX_DIVS.xDiv VOX_DIVS.yDiv]
        
        %RowArr(:,1) = size(crownSegmentsAll,1) - RowArr(:,1);
        %RowArr(:,2) = size(crownSegmentsAll,2) - RowArr(:,2);
        
        if(and(NC.ISPLOTON,plot3D))
            fprintf('\n');
            disp(strcat({'   Printing 3D Features...'}));         
            f2=figure('name','3D Texture Features');
            set(f2, 'Position', [970 620 950 500]);
            title('12 GLCM Features (row wise), for each depth (column wise) ');
        end
        
        voxelFeatureWiseIndices = zeros(VOX_DIVS.xDiv, VOX_DIVS.yDiv, size(crownSegmentsAll,2),VOX_DIVS.zDiv);        
        for levl = 1:VOX_DIVS.zDiv                  
            PC1Mat = zeros(VOX_DIVS.xDiv*VOX_DIVS.yDiv);
            finalFeatureArr = zeros(VOX_DIVS.xDiv*VOX_DIVS.yDiv, VOX_DIVS.zDiv);
            for fCount = 1:1:size(fNames,2)  %numel(fieldnames(fNames))    % fet2skip  % for fCount = 1+fet2skip:size(fNames,1)-1                   
               tempTxtArr = f3Dvoxel.(fNames{fCount})(:,:,levl);  %%%%%
               tempTxtArr = mat2gray(tempTxtArr)*255;
               tempTxtArr(isnan(tempTxtArr))=0;
               % interpolate 
               idx = find(tempTxtArr==0);
               [r,c] = ind2sub(size(tempTxtArr),idx);
               for j = 1:1:size(r,1)
                  avgVal = getRC(r(j),c(j),tempTxtArr);
                  tempTxtArr(r(j),c(j)) = avgVal;
               end               
               tempTxtArr(tempTxtArr<mean(tempTxtArr(:))) = 0;
               
               tempTxtArr = imgaussfilt(tempTxtArr,2);
               
               finalFeatureArr(:,fCount) = tempTxtArr(:);
            end

            COEFF = pca(finalFeatureArr); 
            Itransformed = finalFeatureArr*COEFF;
            PC1Mat = reshape(Itransformed(:,1)+Itransformed(:,2)+Itransformed(:,3),VOX_DIVS.xDiv, VOX_DIVS.yDiv);                    
            
            PC1Mat(1:2,:) = nan; PC1Mat(:,1:2) = nan;
            PC1Mat(end-2:end,:) = nan; PC1Mat(:,end-2:end) = nan;                    

            plotOn = and(NC.ISPLOTON,plot3D);            
            if(plotOn);
                subplot(round(VOX_DIVS.zDiv/4+1,0),4,levl);
            end
            
            PC1Mat(isnan(PC1Mat)) = 0;
            PC1Mat = flipud(PC1Mat);
            
            PC1Mat(isnan(PC1Mat)) = 0;
            PC1Mat = mat2gray(PC1Mat);

            PC1Mat = imresize(PC1Mat, [size(crownSegmentsAll{1},1), size(crownSegmentsAll{1},2)], 'nearest');

            %p
            
            segmentedCSM = treeWatershedCrownSM(PC1Mat*256, ...
            'markers', RowArr, ... % RowArr
            'minHeight', 1, ...
            'fig', false, ...
            'verbose', true);
        
            
           % segmentedCSM = watershedNewNewTestIn(1 - PC1Mat, true); % watershedNew(1 - mat2gray(PC1Mat), seeds);
            voxelFeatureWiseIndices(:,:,:,levl) =  getOverlappedSegment(crownSegmentsAll, segmentedCSM, domTreeThrshLevel, levl,VOX_DIVS);
           % voxelFeatureWiseIndices(:,:,levl,) = segResultd;
           % voxelFeatureWiseIndices(:,:,levl) = imgaussfilt(voxelFeatureWiseIndices(:,:,levl),10);
            
            if(and(NC.ISPLOTON,plot3D))
                imagesc(segmentedCSM); %imagesc(voxelFeatureWiseIndices(:,:,1,levl));
                hold on;
                
                for iCnt = 1:1:size(crownSegmentsAll,2)
                    tmp = crownSegmentsAll{iCnt};
                    tmp(tmp>1) = 1;
                    tmp(:,1) = 0; tmp(:,end) = 0;
                    ImgVector = reshape(edge(tmp), 1, []);
                    [indx,indy] = ind2sub(size(segmentedCSM),find(ImgVector==1));
                    %plot(indy,indx,'.k');
                    plot(indy,indx,'.','MarkerSize',12, 'Color',colorArr(iCnt,:));
                    alpha(0.5);
                end
                
                colormap('jet');
            end
        end
        
end


            %PC1Mat = mat2gray(imgaussfilt(PC1Mat,1));
            %PC1Mat = imgaussfilt(((PC1Mat)),[1 1]);
            
%            se = strel('disk', 3);
%            Ie = imerode(PC1Mat, se);
%            Iobr = imreconstruct(Ie, PC1Mat);
%            Iobrd = imdilate(Iobr, se);
%            Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
%            Iobrcbr = imcomplement(Iobrcbr);
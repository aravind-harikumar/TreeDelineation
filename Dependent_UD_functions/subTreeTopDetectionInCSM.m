function [peaksSubmerged, peakXYOrig, seedImageArr, colorArr, upperEnvelopeXVal, upperEnvelopeYVal, x, y] = ...
    subTreeTopDetectionInCSM(segmentedcsmcdm,nDSM, thresholdBDEM, groundClearencePerc, peakParam, envSmooth, plotOn, ifEllipse)
    shiftVal = 5;
    nDSM = [nDSM(:,end-(shiftVal-1):end) nDSM nDSM(:,1:shiftVal)];
    segmentedcsmcdm = [segmentedcsmcdm(:,end-(shiftVal-1):end) segmentedcsmcdm segmentedcsmcdm(:,1:shiftVal)];

    peaksSubmerged = 0;

    nDSM = mat2gray(flip(nDSM));
    nDSM(:,1) = 0; nDSM(:,end) = 0;
    
    %segmentedcsmcdm = mat2gray(flip(segmentedcsmcdm));
    %segmentedcsmcdm(:,1) = 0; segmentedcsmcdm(:,end) = 0;
    
    % blacken a buffer section (i.e., N rows at the bottom) of the nDSM
    %groundClearence = 1+size(nDSM,1)*groundClearencePerc;
    I_ORG_SIZE = zeros(size(nDSM));
    %nDSM = nDSM(groundClearence:end,:);
    
    % perform otsu thrsholding
    [counts,~] = imhist(nDSM,256);
    thresholdBDEM = otsuthresh(counts);
    %binaryDeM = imbinarize(nDSM, thresholdBDEM); %
    %binaryDeM = imbinarize(imgaussfilt(nDSM,2), 0.1); 
    
    %binaryDeM = binaryDeM(groundClearence:end,:);
    
    % get edge map
    binaryDeM = segmentedcsmcdm;    
    groundClearence = 1+size(binaryDeM,1)*groundClearencePerc;
    edgemap = edge(binaryDeM);
    edgemap(1,:) = 1;
    I_ORG_SIZE(groundClearence:end,:) = edgemap;
    
    % get upper envelope from edgemap
    for ii = 1:1:size(I_ORG_SIZE,2)        
        if(sum(I_ORG_SIZE(:,ii)) > 1)
            [envX, ~] = ind2sub( size(I_ORG_SIZE), find(I_ORG_SIZE(:,ii)==1) );       
            [ff,~] = max(envX);
            I_ORG_SIZE(:,ii) = 0; 
            I_ORG_SIZE(ff,ii) = 1; 
        end
    end
    if(plotOn)
       figure;
       imagesc(binaryDeM);hold on; imagesc(edgemap);
    end
    
    % get envolope indices
    indicex = find(I_ORG_SIZE(:)==1);
    [y,x] = ind2sub( size(I_ORG_SIZE), indicex) ;
    if(plotOn)
        hold on;
        %imagesc(I_ORG_SIZE);
        plot(x,y,'.','Color',[1 0.8 0],'MarkerSize',50); %alpha(0.5);
    end
    
    % smoothen the boundary
    upperEnvelopeXVal = 0:.01:max(x);
    upperEnvelopeYVal = spline(x,y,upperEnvelopeXVal);
    upperEnvelopeYVal = movmean(upperEnvelopeYVal,envSmooth);
    if(plotOn)
        plot(upperEnvelopeXVal,upperEnvelopeYVal,'-','Color',[0 1 0],'LineWidth', 6 );
    end

    [pks,peakIndex] = findpeaks(upperEnvelopeYVal);
    
    %peakIndex = peakfinder(upperEnvelopeYVal,peakParam,0,1, false);
    peakXYOrig = round([upperEnvelopeXVal(peakIndex)' upperEnvelopeYVal(peakIndex)'],0);
    % avoid all max vclues less than 2 percent.
    peakXYOrig = peakXYOrig( peakXYOrig(:,2)>(0.02*size(I_ORG_SIZE,1)) ,:);    
    peakXYOrig = removeRepetitions(peakXYOrig);
    
    if(plotOn)
        hold on;
        plot(peakXYOrig(:,1),peakXYOrig(:,2),'.r', 'MarkerSize', 50);
        stem(peakXYOrig(:,1),peakXYOrig(:,2));

    end
    
    valleyParam = 0;%(max(upperEnvelopeYVal)-min(upperEnvelopeYVal))/12;
    ValleyIndex = peakfinder(upperEnvelopeYVal,valleyParam,100,-1, true);
    ValleyXYOrig = round([upperEnvelopeXVal(ValleyIndex)' upperEnvelopeYVal(ValleyIndex)'],0);
    ValleyXYOrig(ValleyXYOrig(:,2)<0,2) = 0;
    %ValleyXYOrig = ValleyXYOrig( ValleyXYOrig(:,2)>0.8*size(I_ORG_SIZE,1) ,:); 
    if(plotOn)
        hold on;  
        plot(ValleyXYOrig(:,1),ValleyXYOrig(:,2),'.b', 'MarkerSize', 50);
        stem(ValleyXYOrig(:,1),ValleyXYOrig(:,2));  
    end
    
    pointList =  [peakXYOrig; ValleyXYOrig];
    pointListSorted = sortrows(pointList,1);
    
    %nDSM = nDSM(:,shiftVal:end-(shiftVal+1));
    
    set(gca,'Xdir','reverse');
    set(gca,'Ydir','Normal');
    
    seedImageArr = {};
    colorArr = [];
    icnt = 1;
    for iPeakCount = 1:1:size(peakXYOrig,1)
        seedImage = zeros(size(nDSM,1),size(nDSM,2));
        queryRow = peakXYOrig(iPeakCount,:);
        index_A = find(ismember(pointListSorted,queryRow,'rows'));   
        
        % get the left and the right valley points for a selected peak
        % point;
        leftRow = pointListSorted(index_A-1,:);
        rightRow = pointListSorted(index_A+1,:);
        
        % define ellipe parameters
        r_length = rightRow(1)-leftRow(1);
        r_height = queryRow(2);
        
        if( strcmp(ifEllipse, 'ellipse'))
            a=r_length/2; %+ round(size(nDSM,1)*0.05,0); % horizontal radius %+ buffer(for trial)
            b=r_height/2; %+ round(size(nDSM,1)*2,0); % vertical radius %+ buffer(for trial)            
            % x0,y0 ellipse centre coordinates
            x0=queryRow(1);  %x0=queryRow(1)-shiftVal; 
            y0=queryRow(2)/2;
            t=-pi:0.01:pi;
            x=x0+a*cos(t);
            y=y0+b*sin(t);            
            colorVal = rand(1,3);
            if(plotOn)
                plot(x,y,'-.','LineWidth', 5', 'Color', colorVal)
            end
            colorArr = [colorArr; colorVal];            
            [X_Indx,Y_Indx] = meshgrid(1:1:size(nDSM,2),1:1:size(nDSM,1));        
            [in,on] = inpolygon(X_Indx(:),Y_Indx(:),x,y);
            seedImage(or(on,in))=iPeakCount;
            if(plotOn)
                imagesc(seedImage); alpha(0.2);
            end
        else
            rectangle('Position',[leftRow(1) 0 r_length r_height], 'FaceColor', [rand(1) rand(1) rand(1) 0.5]);        
            seedImage(1:queryRow(2),leftRow(1)+1:rightRow(1)) = iPeakCount;
        end     
            tempImg1 = mat2gray(seedImage);
            tempImg1(tempImg1(:)>0) = 1;
            %tempImage = mat2gray(seedImage);
            %tempImage(tempImage<0.5) = 0;
            %nDSM(nDSM>0.1) = 1;
            %nDSM(nDSM<=0.1) = 0;
            %opp = tempImage.*nDSM;
            ss = sum(tempImg1(:));
            %if(ss>1000)
             %   if(ss>(sum(sum(tempImg1))*0.45))
                %if( var(opp(:)) > 0.02 )
             %       ddd = [ddd;ss];
                    seedImageArr{icnt} = flipud(seedImage);
                    icnt = icnt+1;
            %end
         %   end
            %end
    end
end


function [points] = pointsinellipse(xyzdata,bbPosDetails)
    Xvals = xyzdata(:,1);
    Yvals = xyzdata(:,2);
    centerX = bbPosDetails(1) + bbPosDetails(3)/2;
    centerY = bbPosDetails(2) + bbPosDetails(4)/2;
    dist = sqrt(((Xvals(:)-centerX).^2)/((bbPosDetails(3)/2)^2)  + ((Yvals(:)-centerY).^2)/((bbPosDetails(4)/2)^2)); % distance calc.
    %dist = sqrt((Xvals(:)-centerX).^2  + (Yvals(:)-centerY).^2);
    % Find points inside circle
    points = xyzdata(dist<1,:);
end

function extremeLiDARDaraArray = removeRepetitions(extremeLiDARDaraArray)
    DistBwExtrmPoints = tril(squareform(pdist(extremeLiDARDaraArray(:,1),'euclidean')));
    avgDist = median(DistBwExtrmPoints(DistBwExtrmPoints(:)>0));
    [rowIdx,~] = ind2sub(size(DistBwExtrmPoints),find(DistBwExtrmPoints < 10 & DistBwExtrmPoints > 0));
    extremeLiDARDaraArray = removerows(extremeLiDARDaraArray,'ind',unique(rowIdx));    
end


% function [peaksSubmerged, peakXYOrig, seedImage, upperEnvelopeXVal, upperEnvelopeYVal, x, y] = ...
%     testCurveFitTreeTopDetection(nDSM, groundClearencePerc, peakParam, envSmooth, plotOn)
%     
%     %subplot(2,2,2);
%     %nDSM = imgaussfilt(nDSM,5);
%     figure(2);
%     hold on;
%     nDSM = mat2gray(flip(nDSM));
%     nDSM(:,1) = 0; nDSM(:,end) = 0;
%     
%     % blacken a buffer section (i.e., N rows at the bottom) of the nDSM
%     groundClearence = 1+size(nDSM,1)*groundClearencePerc;    
%     I_ORG_SIZE = zeros(size(nDSM));
%     nDSM = nDSM(groundClearence:end,:);
%     
%     % perform otsu thrsholding
%     [counts,~] = imhist(nDSM,8);
%     TreshVal = otsuthresh(counts);
%     binaryDeM = imbinarize(nDSM,TreshVal+0.1);
%     if(plotOn)  
%         dd = zeros(size(binaryDeM));
%         dd(1) =1;
%         imagesc(dd);
%         %imagesc(binaryDeM);
%         %colormap(jet);
%         axis([1 size(I_ORG_SIZE,1) 0 size(I_ORG_SIZE,2)])
%         set(gcf, 'Position',get(0,'Screensize'));
%     end
%     
%     % get edge map
%     edgemap = edge(binaryDeM);
%     edgemap(1,:) = 1;
%     I_ORG_SIZE(groundClearence:end,:) = edgemap;
%     if(plotOn)
%         hold on;
%         %imagesc(edgemap);
%        % alpha(0.5);
%     end
%     
%     % get upper envelope from edgemap
%     for ii = 1:1:size(I_ORG_SIZE,2)        
%         if(sum(I_ORG_SIZE(:,ii)) > 1)
%             [envX, ~] = ind2sub( size(I_ORG_SIZE), find(I_ORG_SIZE(:,ii)==1) );       
%             [ff,~] = max(envX);
%             I_ORG_SIZE(:,ii) = 0; 
%             I_ORG_SIZE(ff,ii) = 1; 
%         end
%     end
%     if(plotOn)
%         %imagesc(binaryDeM);
%         hold on;
%         %imagesc(I_ORG_SIZE);
%         %alpha(0.5);
%     end
%     
%     % get envolope indices
%     indicex = find(I_ORG_SIZE(:)==1);
%     [y,x] = ind2sub( size(I_ORG_SIZE), indicex) ;
%     if(plotOn)
%         imagesc(I_ORG_SIZE);
%         %plot(x,y,'.','Color',[1 1 0],'MarkerSize',50); %alpha(0.5);
%         hold on;
%     end
% 
%     
%     % smoothen the boundary
%     upperEnvelopeXVal = 0:.01:max(x);
%     upperEnvelopeYVal = spline(x,y,upperEnvelopeXVal);
%     plot(upperEnvelopeXVal,upperEnvelopeYVal,'-','Color',[0 1 0],'LineWidth', 6 );
%     upperEnvelopeYVal = movmean(upperEnvelopeYVal,envSmooth);
%    % plot(upperEnvelopeXVal,upperEnvelopeYVal,'-','Color',[1 0 0],'LineWidth', 12 );
%    % plot(x,y,'.','Color',[1 1 0],'MarkerSize',50); %alpha(0.5);
%     % Find peak values
%     % FYI - Origin in matlab is left-top : so size(I_ORG_SIZE,1)-y1
%     peakParam = 1;%(max(upperEnvelopeYVal)-min(upperEnvelopeYVal))/12;
%     peakIndex = peakfinder(upperEnvelopeYVal,peakParam,0,1, false);
%     peakXYOrig = round([upperEnvelopeXVal(peakIndex)' upperEnvelopeYVal(peakIndex)'],0); 
%     peakXYOrig = peakXYOrig( peakXYOrig(:,2)>0.2*size(I_ORG_SIZE,1) ,:); 
%     if(plotOn)
%         plot(peakXYOrig(:,1),peakXYOrig(:,2),'.r', 'MarkerSize', 150);
%         stem(peakXYOrig(:,1),peakXYOrig(:,2));
%         hold on;
%     end
%     
%     valleyParam = 0.5;%(max(upperEnvelopeYVal)-min(upperEnvelopeYVal))/12;
%     ValleyIndex = peakfinder(upperEnvelopeYVal,valleyParam,100,-1, false);
%     ValleyXYOrig = round([upperEnvelopeXVal(ValleyIndex)' upperEnvelopeYVal(ValleyIndex)'],0); 
%     peakXYOrig = peakXYOrig( peakXYOrig(:,2)>0.8*size(I_ORG_SIZE,1) ,:); 
%     if(plotOn)
%         plot(ValleyXYOrig(:,1),ValleyXYOrig(:,2),'.b', 'MarkerSize', 150);
%         stem(ValleyXYOrig(:,1),ValleyXYOrig(:,2));
%         hold on;
%     end
%     
%     % submerged peaks (so that it comes to the center of tree
%     %peakSubFactor = 0.5;
%     %peaksSubmerged = round([peakXYOrig(:,1)  (peakXYOrig(:,2) - (peakXYOrig(:,2))*peakSubFactor) + (peakXYOrig(:,2)*groundClearencePerc)],0);
%     %if(plotOn)
%         %plot(peaksSubmerged(:,1),peaksSubmerged(:,2),'.g', 'MarkerSize', 80);    
%     %end
%     %peaksSubmerged = [peaksSubmerged; [25 45]];  % HARD CODED VALUE TO DETECT BACKGROUND
%     % remove peaks > 0.8 percentage
%     %peaksSubmerged = peaksSubmerged( peaksSubmerged(:,2)>0.2*size(I_ORG_SIZE,1) ,:); 
%     %imagesc(binaryDeM);
%     %axis([1 size(I_ORG_SIZE,1) 0 size(I_ORG_SIZE,2)])
% %     colormap('jet'); %alpha(0.5);
% %     cmap_mod = colormap;
% %     cmap_mod(end,:) = [0.8 0 0];
% %     cmap_mod(1,:) = [0 0 0.7];
% %     colormap(cmap_mod);
%     
% %     seedImage = zeros(size(I_ORG_SIZE));    
% %     indx = sub2ind(size(I_ORG_SIZE), peaksSubmerged(:,2),peaksSubmerged(:,1) );
% %     seedImage(indx) = 1;    
% %     se = strel('line', round(size(I_ORG_SIZE,1)*0.05,0), 90); %strel('diamond', 2); %strel('disk', 1) ; STREL('square',W) ; strel('line', 20, 90);
% %     se1 = strel('diamond', round(size(I_ORG_SIZE,1)*0.01,0));
% %     seedImage = imdilate(seedImage, se)|imdilate(seedImage, se1);
% %     seedImage = imdilate(seedImage, se);
% %     seedImage(seedImage>0) = 1;
% %     if(plotOn)
% %         imagesc(seedImage); alpha(0.5);
% %     end
% %     seedImage = flip(seedImage);
%     
% end
% 
% %     p = polyfit(x,y,20);
% %     upperEnvelopeXVal = linspace(0, size(I_ORG_SIZE,2), length(x));
% %     upperEnvelopeYVal = polyval(p,upperEnvelopeXVal);
% %     if(plotOn)
% %         plot(upperEnvelopeXVal,upperEnvelopeYVal);    
% %     end


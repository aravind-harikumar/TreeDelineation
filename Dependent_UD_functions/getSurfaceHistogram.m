function getSurfaceHistogram(currentCloud)
    figure;
    %subplot(1,2,2);
    hold on;
    getPointDensityCover(currentCloud);
end

function getPointDensityCover(currentCloud)
    MMS = getMaxMins(currentCloud);
    xDi = 8; yDi = 8;
    gXStep = abs(MMS.maxX - MMS.minX)/xDi;
    gYStep = abs(MMS.maxY - MMS.minY)/yDi;
    halfXGridSize  = gXStep/2;
    halfYGridSize  = gXStep/2;
    xDimGrid = size(MMS.minX+halfXGridSize : gXStep : MMS.maxX-halfXGridSize,2);
    yDimGrid = size(MMS.minY+halfYGridSize : gYStep : MMS.maxY-halfYGridSize,2);
    
    levels = fliplr(linspace(MMS.minZ,MMS.maxZ,10));
    csm = zeros(xDimGrid,yDimGrid,length(levels));  
    maxVall = max(currentCloud(:,3));
    
    % Generating csm and cdm
    for i = MMS.maxX-halfXGridSize:-gXStep:MMS.minX+halfXGridSize
        for j = MMS.maxY-halfYGridSize:-gYStep:MMS.minY+halfYGridSize
            cond1 = and(currentCloud(:,1)>i-halfXGridSize, currentCloud(:,1)<=i+halfXGridSize);
            cond2 = and(currentCloud(:,2)>j-halfYGridSize, currentCloud(:,2)<=j+halfYGridSize);
            indGridPoints = find(and(and(cond1,cond2),true));
            currentCloud(indGridPoints,3) = currentCloud(indGridPoints,3)+(maxVall-max(currentCloud(indGridPoints,3)));
        end
    end
            
%     plot3(currentCloud(:,1),currentCloud(:,2),currentCloud(:,3),'*','MarkerSize',10, 'Color', [1 0 0]);
%     axis([min(currentCloud(:,1)) max(currentCloud(:,1)) min(currentCloud(:,2)) max(currentCloud(:,2)) min(currentCloud(:,3)) max(currentCloud(:,3))]);
%     xlabel('Tree Height'); ylabel('Distance to the refernce cyclinder surface'); zlabel('Distance to Stem');

    maxValcsm = 0;
    for levelCnt = 1:1:length(levels)-1
        cnt1 = 1;
        for i = MMS.maxX-halfXGridSize:-gXStep:MMS.minX+halfXGridSize
            cnt2 = 1;
            for j = MMS.maxY-halfYGridSize:-gYStep:MMS.minY+halfYGridSize              
                cond1 = and(currentCloud(:,1)>i-halfXGridSize, currentCloud(:,1)<=i+halfXGridSize);
                cond2 = and(currentCloud(:,2)>j-halfYGridSize, currentCloud(:,2)<=j+halfYGridSize);
                cond3 = and(currentCloud(:,3)<=levels(levelCnt), currentCloud(:,3)>levels(levelCnt+1));
                indGridPoints = find(and(and(cond1,cond2),cond3));
                if(size(indGridPoints,1)>0) 
                    plot3(currentCloud(indGridPoints,1),currentCloud(indGridPoints,2),currentCloud(indGridPoints,3),'*','MarkerSize',10, 'Color', [1 0 0]);
                    axis([min(currentCloud(:,1)) max(currentCloud(:,1)) min(currentCloud(:,2)) max(currentCloud(:,2)) min(currentCloud(:,3)) max(currentCloud(:,3))]);
                    xlabel('Tree Height'); ylabel('Distance to the refernce cyclinder surface'); zlabel('Distance to Stem');
                    camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-90,90); box on;
                    csm(cnt1,cnt2,levelCnt) = numel(currentCloud(indGridPoints,3));                    
                    if(maxValcsm<csm(cnt1,cnt2,levelCnt))
                        maxValcsm = csm(cnt1,cnt2,levelCnt);                    
                    end                    
                end
                cnt2 = cnt2+1;
            end
            cnt1 = cnt1+1;
        end
    end
    
    ii = 1;
    cnt1 = 1; 
    AA = [];
    for i = MMS.maxX-halfXGridSize:-gXStep:MMS.minX+halfXGridSize
        cnt2 = 1;
        for j = MMS.maxY-halfYGridSize:-gYStep:MMS.minY+halfYGridSize
            opt = csm(cnt1,cnt2,:);
            if(and(size(opt(:),1)>0,sum(opt(:))>(maxValcsm*0.25)))
                subplot(xDi,yDi,ii); ii = ii+1;
                plot(1:1:length(levels),normalize(opt(:))','-*');  
                
                sds =[1:1:length(levels)]; vsdf = normalize(opt(:))';
                AA = [AA sds*vsdf'];
                axis([1 length(levels) 0 1]);
                cnt2 = cnt2+1;
            end
        end
        cnt1 = cnt1+1;
    end
    title(num2str(var(AA)))
  
end

function Norm=normalize(A)
    Norm_A=[];
    for i=1:size(A,2)
        Norm_col=(A(:,i)-min(A(:,i)))/(max(A(:,i))-min(A(:,i)));
        Norm_A=[Norm_A,Norm_col];
    end
    Norm=Norm_A;
end
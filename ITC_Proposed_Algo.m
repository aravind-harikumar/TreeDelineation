function [treeTopsArray, X_ExtendArr, Y_ExtendArr, stData] = ITC_Proposed_Algo(fullFileName,outLASFilepathLS,treeCnt)
% close all existing open figures
% close all;
%---------------------------------------------------------------%
%                       Read LiDAR data
%----------------------------------------------------------------%
fprintf('Loading candidate segment data...');
singleTreeLiDARdata = LoadData(fullFileName);
fprintf('done!\n');
% LiDAR data Pre-processing and basic tree information retrieval
fprintf('Preprocessing candidate segment data...');
stData = dataPreProcessing(singleTreeLiDARdata);
fprintf('done!\n');

%----------------------------------------------------------------%
%               Plot LiDAR data (original space)
%----------------------------------------------------------------%
plotStem = true; %FYI: NC.PLOT also should be true for plotting.
plotLiDARData(stData, plotStem);

%----------------------------------------------------------------%
%                 Generate projected LiDAR data
%----------------------------------------------------------------%
fprintf('Projecting candidate segment data...');
scaleFactor = 1; % Fixed
[prjData, DiffDist, EV_ratio_thresh, EVals] = getProjectedLiDARdata(stData, scaleFactor);
saveProjectedDataAsLas(outLASFilepathLS,'out.las',prjData.lidarDataArray);
fprintf('done!\n');
voxDimCalculation = false;

if(EV_ratio_thresh<0.95)
    %----------------------------------------------------------------%%%%
    %                     Get optimal voxel size
    %----------------------------------------------------------------%%%%

    if(voxDimCalculation)
    fprintf('Computing optimal voxel dimension...');
    filterOrder = 2 ; plotCSM = true; %(for testing)
    %VOXEL_DIM_VARIOGRAM.xDiv = 200; VOXEL_DIM_VARIOGRAM.yDiv = 200;
    VOXEL_DIM_FINE = 0.25;
    VOXEL_XYZ_DIVS = getVoxelDimension(stData, VOXEL_DIM_FINE);
    [csmcdmNormalized, dsmImageSegmentedArr, colorArr] = getCSMCDM(EV_ratio_thresh, ...  EV_ratio_thresh
    prjData.lidarDataArray, VOXEL_XYZ_DIVS, filterOrder, plotCSM);
    % Sill value from semi-variagram
    sillVal = getSill(csmcdmNormalized, true)/2;
    fprintf('done!\n'); 
    end

    %----------------------------------------------------------------%%%%
    %                 2D CSM feature generation
    %----------------------------------------------------------------%%%%
    fprintf('Generating CSM and generating approximate crown outline...');
    % Preprocesiing
    OPTIMAL_VOXEL_DIM = 0.1; %sillVal;
    VOXEL_XYZ_DIVS = getVoxelDimension(prjData, OPTIMAL_VOXEL_DIM);
    % Generate dsm image
    filterOrder = 2; plotCSM = true; %(for testing)
    [csmcdmNormalized, dsmImageSegmentedArr, colorArr] = getCSMCDM(EV_ratio_thresh, ...  EV_ratio_thresh
        prjData.lidarDataArray, VOXEL_XYZ_DIVS, filterOrder, plotCSM);
    fprintf('done!\n');
    if(true)
    %----------------------------------------------------------------%%%%
    %                 3D GLCM feature generation
    %----------------------------------------------------------------%%%%
    fprintf('Generating voxel attribute...');
    % generate voxel level attributes
    PRINT_VOXELS = false;
    [f3D_Perf, pdDataParmas] = formVoxelSpacefromTLSProcessing(prjData,...
        VOXEL_XYZ_DIVS,PRINT_VOXELS);
    pdDataParmas.lidarDataArrayComplete = pdDataParmas.lidarDataArray;
    fprintf('done!\n');
    
    fprintf('Generating 3D GLCM features...');
    % Get voxel GLCM features
    plot3D = true;
    % get threshold to obtain dominant tree in the projected space
    [counts,~] = imhist(csmcdmNormalized,256);
    thresholdfromcsmcdm = otsuthresh(counts);
    % get voxel GLCM features
    bias = 0.25;
    domTreeThrshLevel = VOXEL_XYZ_DIVS.zDiv*(DiffDist+bias); %(1-DiffDist); %%%%%%%%%%%%%%%%%%% % test % DiffDist+0.2
    [~, subTreeVoxelsByFeature] =  getvoxelGLCMFeatures(f3D_Perf, VOXEL_XYZ_DIVS, plot3D, ...
        dsmImageSegmentedArr, domTreeThrshLevel, colorArr);
    fprintf('done!\n');
    
    %----------------------------------------------------------------%
    %               Plot 3D voxel plot (projected space)
    %----------------------------------------------------------------%
    printVoxel = true;
    colorArr = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0; 0.5 0.5 0; 0 0.5 0; 0 0 0.3; 0 0.3 0; 0 0.5 0.4];
    if(printVoxel)
        fprintf('Printing voxel space...');
    end
    formVoxelSpace(subTreeVoxelsByFeature, pdDataParmas, VOXEL_XYZ_DIVS, printVoxel,colorArr);
    if(printVoxel)
        fprintf('done!\n');
    end
    %----------------------------------------------------------------%%%%
    %          SubDominant and Dominant Tree Deta Extraction
    %----------------------------------------------------------------%%%%
    fprintf('Detecting crown top and extend...');
    plotResult = true;
    [treeTopsArray,X_ExtendArr,Y_ExtendArr] = plotCluteredTrees(stData, pdDataParmas, ...
        subTreeVoxelsByFeature, f3D_Perf.voxContainedPointIndice, plotResult,colorArr);
        fprintf('done!\n');
    %crownDIm = max(X_ExtendArr,Y_ExtendArr);
    end
    treeTopsArray=[0 0]; X_ExtendArr=[0 0]; Y_ExtendArr=[0 0];
else
    treeTopsArray=[0 0]; X_ExtendArr=[0 0]; Y_ExtendArr=[0 0];
    fprintf('No subdominant treees detected...');
end

end
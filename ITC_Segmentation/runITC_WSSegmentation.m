% you can get details of the algorirhm at:
% http://mparkan.github.io/Digital-Forestry-Toolbox/tutorial-2.html

% Data Prerequisites:
% The data should contain only 2 classes (i.e., ground, and non ground). If
% the data contain more classes,Run Generate2ClassData.m to generate two class data.
% the terrain class is assumed to be of value 2, and the other class as 1.

function peaks_crh = runITC_WSSegmentation(inputLASPlotPath, outLASFilePath, cutTrees)

    % set paths and variables
    files = dir( fullfile( inputLASPlotPath,'twoclass_000416_elli_elli_Tile_14.las' ));
    for file = files'
        close all;
        inParams.InputFolder = inputLASPlotPath;
        inParams.OutputFolder = outLASFilePath;
        inParams.PlotFileName = file.name;
        inParams.FileType = '';
        % params for the ITC_Algo
        inParams.cellResolution = 0.2;
        inParams.PlotRadius = 111; % should be same as the original plot radius;
        inParams.TreetopDetectionAlgo = 'allometricRadius'; %'allometricRadius','LevelSet','hMaxima','fixedRadius'
        inParams.cutTrees = cutTrees;
        % plot on/off parameters
        inParams.Plot_CHM = false; inParams.Plot_Peaks = false;
        inParams.Plot_WSeG = true; inParams.Plot_3DCluster = true;
        inParams.writeDatatoXLS = true;
        % reference data
%         inParams.ShpFileFolder = ShpFileFolder;
        % run tree detection algorithm
        peaks_crh = Run_ITC_Algo(inParams);
        %figure; plot(peaks_crh(peaks_crh(:,4)==1,1), peaks_crh(peaks_crh(:,4)==1,2), '*'); title('Estimated Crown Positions');
        pause(1);
    end
   
end

function peaks_crh = Run_ITC_Algo(inParams)

pc = LASread( fullfile(inParams.InputFolder, strcat(inParams.PlotFileName,inParams.FileType)) );
xyz = [pc.record.x, pc.record.y, pc.record.z];
classification = pc.record.classification;

[models, refmat] = elevationModels(xyz, ...
    classification, ...
    'classTerrain', [2], ...
    'classSurface', unique(classification)', ...
    'cellResolution', inParams.cellResolution, ...
    'closing', 10, ...
    'smoothingFilter', fspecial('gaussian', [11,11], 1), ...
    'outputModels', {'terrain', 'surface', 'height'}, ...
    'fig', inParams.Plot_CHM, ...
    'verbose', true);

[peaks_crh, ~] = canopyPeaks(models.height.values, ...
    refmat, ...
    'method', inParams.TreetopDetectionAlgo, ... % 'allometricRadius'
    'allometry', @(h) 0.5 + 0.25*log(max(h,1)), ...
    'adjacency', @(h) min(0.5 + 0.5*log(max(h,1)),4), ...
    'minPeakHeight', 255*15/100, ... % height threshold
    'plotRadius', inParams.PlotRadius',...
    'cellresolution', inParams.cellResolution',...
    'fig', inParams.Plot_Peaks, ...
    'verbose', true);

borderMargin = 5;
se = strel('disk', ...
    ceil(borderMargin/(inParams.cellResolution*10)), ...
    6);
mask = imerode(padarray(models.mask, [1 1], false), se);
mask = mask(2:end-1,2:end-1);

label_2d = treeWatershed(models.height.values, ...
    'markers', peaks_crh, ...
    'minHeight', 255*24/100, ... % height threshold
    'mask', mask, ...
    'fig', inParams.Plot_WSeG, ...
    'plotRadius', inParams.PlotRadius',...
    'cellresolution', inParams.cellResolution',...
    'verbose', true);

% convert image to map coordinates, [x y] = [row col 1] * R
[nrows, ncols] = size(label_2d);
[col_grid, row_grid] = meshgrid(1:ncols, 1:nrows);
xy_grid = [row_grid(:), col_grid(:), ones(numel(label_2d),1)] * refmat;

% index vegetation points
idxl_veg = ismember(classification, [1,2]);

label_vec = label_2d(:);
idxl_labelled = (label_vec ~=0);
label_vec = label_vec(idxl_labelled);

% compute nearest neighbours
[idxn_knn, d_nn] = knnsearch(xy_grid(idxl_labelled,:), ...
    xyz(idxl_veg, 1:2), ...
    'K', 1);

label_veg = label_vec(idxn_knn);
label_veg(d_nn > inParams.cellResolution) = 0;
label_3d = zeros(size(xyz,1),1);
label_3d(idxl_veg) = label_veg;
[label_3d(label_3d ~= 0), ~] = grp2idx(label_3d(label_3d ~= 0));

[color, rgb, cmap] = clusterColor(xyz, ...
    label_3d, ...
    'adjacency', '2d', ...
    'buffer', 10, ...
    'colormap', 'cmap12', ...
    'unlabelledColor', [0.1, 0.1, 0.1], ...
    'fig', inParams.Plot_3DCluster, ...
    'verbose', true);

figure; hold on;
lbls = unique(label_3d);
maxTip = [];
for i = 2:1:size(lbls,1)
    clusArr = xyz(label_3d==lbls(i),:);
    [maxval,idx] = max(clusArr(:,3));
   % [minval,idxmin] = min(clusArr(:,3));
     maxTip = [maxTip; clusArr(idx,1:2)];
%      g = fspecial('gaussian', 200, 0.25)*(maxval-minval);
%      surf(g);hold on;
end
plot(maxTip(:,1), maxTip(:,2),'.b', 'MarkerSize', 20); 
title('Estimated Crown Positions');

bbox = [min(xyz(:,1)) min(xyz(:,2)); max(xyz(:,1)) max(xyz(:,2))];


% Save Individual tree data, in the plot, to the outFolder;
% SaveTreedata(inParams,label_3d,label_2d,pc,xyz,peaks_crh);

end

function SaveTreedata(inParams,label_3d,label_2d,pc,xyz,peaks_crh)
% read reference data
bbox = [min(xyz(:,1)) min(xyz(:,2)); max(xyz(:,1)) max(xyz(:,2))];
S = shaperead( fullfile(inParams.ShpFileFolder,'ALL'),'BoundingBox',bbox); % 2013_Alberi_PNA_SerFo
S1 = shaperead( fullfile(inParams.ShpFileFolder,'ALcomm'),'BoundingBox',bbox); % 2013_Alberi_PNA_SerFo
S2 = shaperead( fullfile(inParams.ShpFileFolder,'ALOmm'),'BoundingBox',bbox); % 2013_Alberi_PNA_SerFo

%shapewrite(table2struct(T), '26995_12710_seg_metrics.shp');
%mapshow(S);

% Save data to file
STable = struct2table(S);
STable1 = struct2table(S1);
STable2 = struct2table(S2);
hold on;
plot(table2array(STable(:,2)),table2array(STable(:,3)),'.r', 'MarkerSize', 20);
hold on;
plot(table2array(STable1(:,2)),table2array(STable1(:,3)),'.g', 'MarkerSize', 20);
hold on;
plot(table2array(STable2(:,2)),table2array(STable2(:,3)),'.b', 'MarkerSize', 20);

%figure;hist(table2array(STable(:,8)));
if(inParams.cutTrees)
% Save individual tree point cloud data
treeArrTable = {};
IndividualTreeID = unique(label_3d);
for i = 2:1:length(IndividualTreeID)
    idxl_sample = (label_3d == IndividualTreeID(i));
    [getTreeAttr,retMaxXYZ,outFileName] = GetTreeData(pc, inParams.PlotFileName, IndividualTreeID(i), idxl_sample, inParams.OutputFolder);
    treeArrTable{i,1} = outFileName;
    treeArrTable{i,2} = getTreeAttr.header.max_z;
    
    % get closest reference point
    DetectedPeak = [retMaxXYZ(1) retMaxXYZ(2)];
    refTreeLocationsAll = table2array([STable(:,2) STable(:,3)]);
    % find nearest reference tree XY (else print [0 0])
    [Idx,D] = knnsearch(refTreeLocationsAll,DetectedPeak);
    if(D<2)
        treeArrTable{i,3} = refTreeLocationsAll(Idx,1);
        treeArrTable{i,4} = refTreeLocationsAll(Idx,2);
    else
        treeArrTable{i,3} = 0;
        treeArrTable{i,4} = 0;
    end
    
    % detected tree top location
    treeArrTable{i,5} = DetectedPeak(1);
    treeArrTable{i,6} = DetectedPeak(2);
    
    % rmse in tree location
    treeArrTable{i,7} = pdist([refTreeLocationsAll(Idx,1:2);DetectedPeak]);
    treeArrTable{i,8} = str2double(cell2mat(table2array(STable(Idx,13))));
    treeArrTable{i,9} = getMaxCrownWidthFromSegment(IndividualTreeID(i), label_2d, peaks_crh, inParams.PlotRadius);
    % get crown area
end

% save data to xls
if(inParams.writeDatatoXLS)
    treeArrTable = cell2table(treeArrTable);
    % label the columns
    treeArrTable.Properties.VariableNames([1 2 3 4 5 6 7 8 9]) = {'FileName' 'Height' 'X_REFERENCE' 'Y_REFERENCE' 'X_DETECTED' 'Y_DETECTED' 'RMSE' 'DBH' 'MAXRADIUS'};
    WriteData2XLS(treeArrTable,'DomTree.xls');
end

end
end

function scaledMaxCrownRadius = getMaxCrownWidthFromSegment(indx, label_3d, peaks_crh, PlotRadius)

idxl_sample = (label_3d == indx);

area = regionprops(idxl_sample , 'Area');
if( and(~isempty(area),size(area,1)==1) )
    if(area.Area>10)
        segboundary = find(edge(idxl_sample));
        [Boundary_X, Boundary_Y] = ind2sub(size(label_3d),segboundary);
        in = inpolygon(peaks_crh(:,2),peaks_crh(:,1),Boundary_X,Boundary_Y);
        BoundaryPoints = [Boundary_X, Boundary_Y];
        ContainedPeak = [peaks_crh(in,2) peaks_crh(in,1)];

        % dist b/w peak and boundary points
        if(and(~isempty(ContainedPeak), ~isempty(BoundaryPoints)))
            D = pdist2(BoundaryPoints,ContainedPeak(1,:));
            [MaxCrownRadius , ~] = max(D);
            scaledMaxCrownRadius = (PlotRadius/256)*MaxCrownRadius;
        else
            scaledMaxCrownRadius = 0;
        end
    end
else
    scaledMaxCrownRadius = 0;
end
%figure; plot(Boundary_X,Boundary_Y,'*');
%axis([-255 255 -255 255])

end


function [r,retMaxXYZ,outFileName] = GetTreeData(pc, PlotFileName, treeID, idxl_sample, outputFolder)

% set record format
switch pc.header.point_data_format_id
    case 1 % 1 -> 3
        recordFormat = 3;
    case 4 % 4 -> 5
        recordFormat = 5;
    case 6 % 6 -> 7
        recordFormat = 7;
    case 9 % 9 -> 10
        recordFormat = 10;
    otherwise % 2,3,5,7,8,10
        recordFormat = pc.header.point_data_format_id;
end

% get the record as it is
r = pc;
r.record.x = pc.record.x(idxl_sample);
r.record.y = pc.record.y(idxl_sample);
r.record.z = pc.record.z(idxl_sample);
r.record.intensity = pc.record.intensity(idxl_sample);
r.record.return_number = pc.record.return_number(idxl_sample);
r.record.number_of_returns = pc.record.number_of_returns(idxl_sample);
r.record.scan_direction_flag = pc.record.scan_direction_flag(idxl_sample);
r.record.flightline_edge_flag = pc.record.flightline_edge_flag(idxl_sample);
r.record.classification = pc.record.classification(idxl_sample);
r.record.classification_synthetic = pc.record.classification_synthetic(idxl_sample);
r.record.classification_keypoint = pc.record.classification_keypoint(idxl_sample);
r.record.classification_withheld = pc.record.classification_withheld(idxl_sample);
r.record.scan_angle = pc.record.scan_angle(idxl_sample);
r.record.user_data = pc.record.user_data(idxl_sample);
r.record.point_source_id = pc.record.point_source_id(idxl_sample);
r.record.gps_time = pc.record.gps_time(idxl_sample);
try
r.record.red = pc.record.red(idxl_sample);    
r.record.green = pc.record.green(idxl_sample);
r.record.blue = pc.record.blue(idxl_sample);
catch
    printf('Error while RGB values! \n')
end
% change header details
r.header.n_point_records = sum(idxl_sample);
r.header.n_point_records_extended = sum(idxl_sample);
r.header.max_x = max(pc.record.x(idxl_sample));
r.header.min_x = min(pc.record.x(idxl_sample));
r.header.max_y = max(pc.record.y(idxl_sample));
r.header.min_y = min(pc.record.y(idxl_sample));
r.header.max_z = max(pc.record.z(idxl_sample));
r.header.min_z = min(pc.record.z(idxl_sample));


lidarDataAr = [r.record.x r.record.y r.record.z];
[retMaxXYZ, ~] = findMaxHeightXY(lidarDataAr);

[~,NAME,~]  = fileparts(PlotFileName);
outFileName = strcat(NAME,'_tree',num2str(treeID),'.las');
LASwrite(r, ...
    fullfile(outputFolder, outFileName), ...
    'version', 13, ...
    'systemID', 'SEGMENTATION', ...
    'recordFormat', recordFormat, ...
    'verbose', true);

end
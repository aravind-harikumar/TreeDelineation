function Generate2ClassData()
    InputFolder = 'C:\My_Files\1_PhD_Research\1_PhD_Research_Topic\1_LiDAR_Forest_Applications\2_Matlab_Tools_LiDAR_Tree_Delination\1_LiDARTileIndexByLocation\Cut_Plots';
    OutputFolder = '.';
    allFiles = dir(fullfile(InputFolder,strcat('*.las')));
    for plotFile = allFiles' 
        GeneratetwoClassData(InputFolder,OutputFolder,plotFile.name);
    end
end

function GeneratetwoClassData(InputFolder,OutputFolder,datasetname)
pc = LASread( fullfile(InputFolder, datasetname) );

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
%cIndx = classification~=2;
newClass = pc.record.classification(:);
newClass(newClass~=2) = 1;
r.record.x = pc.record.x(:);
r.record.y = pc.record.y(:);
r.record.z = pc.record.z(:);
r.record.intensity = pc.record.intensity(:);
r.record.return_number = pc.record.return_number(:);
r.record.number_of_returns = pc.record.number_of_returns(:);
r.record.scan_direction_flag = pc.record.scan_direction_flag(:);
r.record.flightline_edge_flag = pc.record.flightline_edge_flag(:);
r.record.classification = newClass; %pc.record.classification(:);
r.record.classification_synthetic = pc.record.classification_synthetic(:);
r.record.classification_keypoint = pc.record.classification_keypoint(:);
r.record.classification_withheld = pc.record.classification_withheld(:);
r.record.scan_angle = pc.record.scan_angle(:);
r.record.user_data = pc.record.user_data(:);
r.record.point_source_id = pc.record.point_source_id(:);
r.record.gps_time = pc.record.gps_time(:);
%r.record.red = pc.record.red(:);    
%r.record.green = pc.record.green(:);
%r.record.blue = pc.record.blue(:);

% change header details
%r.header.n_point_records = sum(cIndx);

LASwrite(r, ...
    fullfile(OutputFolder, strcat('twoclass_',datasetname) ), ...
    'version', 13, ...
    'systemID', 'SEGMENTATION', ...
    'recordFormat', recordFormat, ...
    'verbose', true);

end
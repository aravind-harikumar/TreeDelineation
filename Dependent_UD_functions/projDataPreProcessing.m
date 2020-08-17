function dataAtt = projDataPreProcessing(inputData)
    %inputData = flattenProjectedCloud(inputData);

    dataAtt.lidarDataArray = inputData;
    dataAtt.lidarDataArrayComplete = inputData;

    dataAtt.XMax = max(inputData(:,1));
    dataAtt.XMin = min(inputData(:,1));

    dataAtt.YMax = max(inputData(:,2));
    dataAtt.YMin = min(inputData(:,2));

    dataAtt.ZMax = max(inputData(:,3));
    dataAtt.ZMin = min(inputData(:,3));
end
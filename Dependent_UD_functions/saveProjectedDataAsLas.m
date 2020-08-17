function saveProjectedDataAsLas(outFilepath,FileName,slabbedCoordinates)
    %dataAtt = write2tablenew(slabbedCoordinates);
    writetable(array2table(slabbedCoordinates), fullfile(outFilepath, 'TempFile.txt'))
    exePath = 'C:\My_Files\3_Software_Packages\Lidar_Tools\LiDAR_CPP_Tools\LasTools\bin\';
    isSuccess = makeLASFromTXT(exePath,'txt2las', strcat(outFilepath,'\'),'TempFile.txt', strcat(outFilepath,'\'),FileName);
end
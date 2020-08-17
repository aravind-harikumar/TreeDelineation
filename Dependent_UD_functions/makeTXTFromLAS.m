function makeTXTFromLAS(lidarDataArray,fullFileName,format)
    dlmwrite(strcat(fullFileName,format),lidarDataArray,'delimiter',' ','precision',10);
end
function makeTXTFromLASMAIN(lidarDataArray,outFilepath)        
    makeTXTFromLAS(lidarDataArray, strcat(outFilepath,'TempFile'),'.txt');
end

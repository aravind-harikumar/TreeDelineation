function isSuccess = makeLASFromTXT(exePath,exeName,inFilePath,inFileName,outFilePath,outFileName)
    command = char(strcat(exePath,exeName,{' -i '},inFilePath,inFileName,{' -o '},outFilePath,outFileName));
    [status,cmdout] = system(command); 
    isSuccess = status;
end
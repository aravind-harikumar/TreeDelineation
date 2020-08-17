function isSuccess = makeDEMFromLAS(exePath,exeName,inFilePath,inFileName,outFilePath,outFileName, demRes)
    command = char(strcat(exePath,exeName,{' -i '},inFilePath,inFileName,{' -o '},outFilePath,outFileName,{' -step 0.25 '})); %-keep_class 1 2 3 4 5 6 7
    [status,cmdout] = system(command); 
    isSuccess = status;
end
function isSuccess = generateDEMfromCSM(slabbedCoordinates, demRes, exePath, outFilepath)
    makeTXTFromLASMAIN(slabbedCoordinates,outFilepath);      
    isSuccess =  makeLASFromTXT(exePath,'txt2las',outFilepath,'TempFile.txt',outFilepath,'Out.las');
    if(isSuccess == 0)        
        fclose('all');
        delete(strcat(outFilepath,'TempFile.txt'));
       % h = msgbox('Operation Completed Successfully','Success');
    end      
    isSuccess = ~makeDEMFromLAS(exePath,'las2dem',outFilepath,'Out.las',outFilepath,'ABC.tif', demRes);
end
function WriteData2XLS(treeArrTable,fileName)
[fid, ~] = fopen(fileName,'a');
if(fid ~= -1)
    fclose(fid);
    delete(fileName);
    writetable(treeArrTable,fileName);
else
    fprintf('Error writing! Please close the Excel file');
end
end
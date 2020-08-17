function getCRplotGraph
    close all; clear all;
    % Get Proposed method tree count;
    csvData = xlsread('C:\Users\harikumar\Desktop\plots1_PropCR.csv');  
    PropMethDBH = csvData(1:23,3);
    PropMethDBHEsti = csvData(1:23,10);

    plot(PropMethDBHEsti,PropMethDBH,'.','MarkerSize', 25)
    
    title('Tree Detection Performance');
    xlabel('DBH (cm)');
    ylabel('Tree count');
    
    
end
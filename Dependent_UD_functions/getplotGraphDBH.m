function getplotGraphDBH
    close all; clear all;
    % Get Proposed method tree count;
    csvData = xlsread('C:\Users\harikumar\Desktop\AA.xlsx');  
    PropMethDBH = csvData(1:197,1);
    PropMethDBHEsti = csvData(1:197,2);
    PropMethDBHEstiSoA = csvData(1:197,3);

    plot(PropMethDBHEsti,PropMethDBH,'.g','MarkerSize', 25)
    hold on;
    plot(PropMethDBHEstiSoA,PropMethDBH,'.b','MarkerSize', 25)
    
    title('Tree Detection Performance');
    xlabel('Reference DBH (cm)');
    ylabel('Estimated DBH (cm)');
    
    hold on;
    aa= [PropMethDBH PropMethDBHEsti];
    X = [ones(length(aa(:,2)),1) aa(:,2)];
    b = X\aa(:,1);
    yCalc2 = X*b;
    plot(aa(:,2),yCalc2,'g', 'LineWidth',1)
    hA = gcf;
    %set(hA,'GridAlpha',0.8); 
    set(hA,'defaultaxesfontname','Times New Roman');
    axis square; grid on;
   % set(gcf, 'Position', get(0, 'Screensize'))
    
    hold on;
    aa= [PropMethDBH PropMethDBHEstiSoA];
    X = [ones(length(aa(:,2)),1) aa(:,2)];
    b = X\aa(:,1);
    yCalc2 = X*b;
    plot(aa(:,2),yCalc2,'b', 'LineWidth',1)
    hA = gcf;
    %set(hA,'GridAlpha',0.8); 
    set(hA,'defaultaxesfontname','Times New Roman');
    axis square; grid on;
    %set(gcf, 'Position', get(0, 'Screensize'))
    set(gca,'fontsize',24)
%     
%     
%     P = polyfit(PropMethDBH,PropMethDBHEsti,1);
%     yfit = P(1)*PropMethDBH+P(2);
%     hold on;
%     plot(PropMethDBH,yfit,'b-.');
%     
%     P = polyfit(PropMethDBH,PropMethDBHEstiSoA,1);
%     yfit = P(1)*PropMethDBH+P(2);
%     hold on;
%     plot(PropMethDBH,yfit,'g-.');
end
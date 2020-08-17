function calcDistanceBasedPerformance()
%SUBDOM TREE DATA
[~,~,raw] = xlsread('C:\Users\harikumar\Desktop\New Foglio di lavoro di Microsoft Excel.xlsx','SUBDOM');  
for idx = 1:numel(raw)
    if isnumeric(raw{idx})
        raw{idx} = num2str(raw{idx});
    end
end
raw(cellfun('isempty',raw))={0};
subTreeData = (raw(1:46,1:8));

%DOM TREE DATA
[~,~,rawD] = xlsread('C:\Users\harikumar\Desktop\New Foglio di lavoro di Microsoft Excel.xlsx','DOM');  
for idx = 1:numel(rawD)
    if isnumeric(rawD{idx})
        rawD{idx} = num2str(rawD{idx});
    end
end
rawD(cellfun('isempty',rawD))={0};
DomTreeData = (rawD(1:163,1:8));

% calcaulte pairwise distance for all trees
op = pdist2(double( str2double(subTreeData(:,3:4)) ),double( str2double(DomTreeData(:,3:4))), 'euclidean');
op(op>25) = 0;
op(op==0) = inf;

TotEstALL = []; TotRealALL = [];
% find tree pairs >=1m and pairs < 2.5m 
[row,col] =  find(and(op>=1,op<2.5));
ReadDBHofSubTrees1mClass = subTreeData(row,2);
EstiDBHofSubTreesw1mClass = subTreeData(row,8);
ReadDBHofDomTrees1mClass = DomTreeData(col,2);
EstiDBHofDomTreesw1mClass = DomTreeData(col,8);
TotReal = str2double([ReadDBHofSubTrees1mClass; ReadDBHofDomTrees1mClass]);
TotEst = str2double([EstiDBHofSubTreesw1mClass; EstiDBHofDomTreesw1mClass]);
notnanRows = ~or(isnan(TotEst), isnan(TotReal));
TotReal = TotReal(notnanRows,:);
TotEst = TotEst(notnanRows,:);
%TotEstALL = [TotEstALL; TotEst];
%TotRealALL = [TotRealALL; TotReal];
AE = sum(TotEst-TotReal)/length(TotReal)
MAE = sum(abs(TotEst-TotReal))/length(TotReal)
RMSE = sqrt(sum((TotEst-TotReal).^2)/length(TotReal))

subplot(1,3,1);
Xla = TotEst-TotReal;
boxplot(Xla,'Color', 'b', 'Symbol', 'kx');grid on;
hold on
plot(mean(Xla), '.b')
axis([0.6 1.2 -20 20]);

% find tree pairs >=2.5m and pairs < 5m 
[row,col] =  find(and(op>=2.5,op<5));
ReadDBHofSubTrees1mClass = subTreeData(row,2);
EstiDBHofSubTreesw1mClass = subTreeData(row,8);
ReadDBHofDomTrees1mClass = DomTreeData(col,2);
EstiDBHofDomTreesw1mClass = DomTreeData(col,8);
TotReal = str2double([ReadDBHofSubTrees1mClass; ReadDBHofDomTrees1mClass]);
TotEst = str2double([EstiDBHofSubTreesw1mClass; EstiDBHofDomTreesw1mClass]);
notnanRows = ~or(isnan(TotEst), isnan(TotReal));
TotReal = TotReal(notnanRows,:);
TotEst = TotEst(notnanRows,:);
%TotEstALL = [TotEstALL; TotEst];
%TotRealALL = [TotRealALL; TotReal];
AE = sum(TotEst-TotReal)/length(TotReal)
MAE = sum(abs(TotEst-TotReal))/length(TotReal)
RMSE = sqrt(sum((TotEst-TotReal).^2)/length(TotReal))

subplot(1,3,2);
Xar = TotEst-TotReal;
boxplot(Xar,'Color', 'b', 'Symbol', 'kx');grid on;
hold on
plot(mean(Xar), '.b')
axis([0.6 1.2 -20 20]);

% find tree pairs >=4m and pairs < 6m 
[row,col] =  find(and(op>=5,op<8));
ReadDBHofSubTrees1mClass = subTreeData(row,2);
EstiDBHofSubTreesw1mClass = subTreeData(row,8);
ReadDBHofDomTrees1mClass = DomTreeData(col,2);
EstiDBHofDomTreesw1mClass = DomTreeData(col,8);
TotReal = str2double([ReadDBHofSubTrees1mClass; ReadDBHofDomTrees1mClass]);
TotEst = str2double([EstiDBHofSubTreesw1mClass; EstiDBHofDomTreesw1mClass]);
notnanRows = ~or(isnan(TotEst), isnan(TotReal));
TotReal = TotReal(notnanRows,:);
TotEst = TotEst(notnanRows,:);
%TotEstALL = [TotEstALL; TotEst];
%TotRealALL = [TotRealALL; TotReal];
AE = sum(TotEst-TotReal)/length(TotReal)
MAE = sum(abs(TotEst-TotReal))/length(TotReal)
RMSE = sqrt(sum((TotEst-TotReal).^2)/length(TotReal))

subplot(1,3,3);
Xab = TotEst-TotReal;
boxplot(Xab, 'Color', 'b', 'Symbol', 'kx');grid on;
axis([0.6 1.2 -20 20]);
hold on
plot(mean(Xab), '.b')
hold off

% %ALL TREE DATA
% [~,~,rawAll] = xlsread('C:\Users\harikumar\Desktop\plots1_PropCR.csv');  
% for idx = 1:numel(rawAll)
%     if isnumeric(rawAll{idx})
%         rawAll{idx} = num2str(rawAll{idx});
%     end
% end
% rawAll(cellfun('isempty',rawAll))={0};
% allData = str2double((rawAll(1:197,3:6)));
% TotReal = allData(:,1);
% TotEst = allData(:,2);
% AE = sum(TotEst-TotReal)/length(TotReal)
% MAE = sum(abs(TotEst-TotReal))/length(TotReal)
% RMSE = sqrt(sum((TotEst-TotReal).^2)/length(TotReal))
% 
% 

end
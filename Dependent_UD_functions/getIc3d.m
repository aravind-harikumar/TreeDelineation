
function retIndicesArr =  getIc3d(lendistArr, lendirecArr, inc)
    retIndicesArr = [];   
    for i=1:1:12    
        retIndices=[];
        for d1Ins = 0:1:lendistArr-1
            d1 = d1Ins*(lendirecArr*inc);            
            for j= 0:1:lendirecArr-1
                retIndices = [retIndices inc*j+(i+d1)];
            end
        end
        retIndicesArr = [retIndicesArr;retIndices];
    end    
end
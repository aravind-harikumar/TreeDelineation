function csm = interploteImage(csm, method)
idx = find(csm==0);
[r,c] = ind2sub(size(csm),idx);
for j = 1:1:size(r,1)
    avgVal = getRC(r(j),c(j),csm);
    csm(r(j),c(j)) = avgVal;
end
if(strcmp(method,'Mean'))
    csm(csm<mean(csm(:))) = 0; %mean(csm(:)) %0 or mean value;
elseif(strcmp(method,'Median'))
    csm(csm<median(csm(:))) = 0; %median(csm(:)) %0 or mean value
else
    fprintf('no thresholding performed')
end
function imgHseed =LevelSet_Seed_Extraction(img_max,img_max_filtered,MinDistance,height_threshold)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

img_seed = imregionalmax(img_max_filtered);
imgHseed=img_seed.*img_max;

img_seed = bwlabel(img_seed);
%% a pixel per seed
for i_seed = unique(nonzeros(img_seed))'
    
    imgHseedTmp = (img_seed==i_seed).*imgHseed;
    imgHseed(imgHseedTmp~=max(max(imgHseedTmp))&(imgHseedTmp~=0))=0;
        if size(find(imgHseed.*(img_seed==i_seed)>0),1)>1 
            pixels = find(imgHseed.*(img_seed==i_seed)>0);
            pixels(pixels==min(pixels))=[];
            imgHseed(pixels)=0;
        end
    
end
imgHseed(imgHseed<height_threshold)=0;
img_seed = img_seed.*(imgHseed>0);
img_seed = bwlabel(img_seed);

%% min Distance among seed
img_seed = imdilate(img_seed,ones(MinDistance,MinDistance));
imgHseed = (img_seed>0).*img_max;

for i_seed = unique(nonzeros(img_seed))'
    
    imgHseedTmp = (img_seed==i_seed).*imgHseed;
    imgHseed(imgHseedTmp~=max(max(imgHseedTmp))&(imgHseedTmp~=0))=0;
        if size(find(imgHseed.*(img_seed==i_seed)>0),1)>1 
            pixels = find(imgHseed.*(img_seed==i_seed)>0);
            pixels(pixels==min(pixels))=[];
            imgHseed(pixels)=0;
        end
    
end
end
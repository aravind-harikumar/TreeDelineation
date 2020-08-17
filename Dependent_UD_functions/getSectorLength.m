function [sectorLength, disttocentre0] = getSectorLength(refPoint, inputPoint, centrePoint, actualCenter, scaleFactor)
    disttocentre0 = pdist2(inputPoint(:,1:2),actualCenter);
    disttocentre = pdist2(refPoint(:,1:2),actualCenter);    
    refAxisPointXY = getPointAtDistance(refPoint(:,1:2),centrePoint,disttocentre);    
    angleBwVextors = getAngleBetweenVectors(refAxisPointXY, inputPoint(:,1:2));
    sectorLength = 2*pi*scaleFactor*(angleBwVextors/360);
end
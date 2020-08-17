function retPoint = getPointAtDistance(point1,point2,distVal)
    d = sqrt((point1(2)-point1(1))^2 + (point2(2) - point1(2))^2);
    r = distVal / d;
    x3 = r * point2(1) + (1 - r) * point1(1); % find point that divides the segment
    y3 = r * point2(2) + (1 - r) * point1(2); % into the ratio (1-r):r
    retPoint = [x3 y3];
end
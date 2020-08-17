function testPoissonDist()
    data = csvread('C:\My_Files\1_PhD_Research\1_PhD_Research_Topic\1_LiDAR_Forest_Applications\3_Research_Files\3_SubdominatTreeDetection\R\Mehtätalo\Plot20.csv');
    [lambdahat,lambdaci] = poissfit(data,0.05);  
    [H,P,STATS] = chi2gof(data,'cdf',@(z)poisscdf(z,lambdahat),'nbins',7)
end
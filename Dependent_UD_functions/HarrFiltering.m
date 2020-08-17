function [A1Coeff,A2Coeff] = HarrFiltering(filtImage)
    [c,s]=wavedec2(filtImage,2,'haar');

    [H1,V1,D1] = detcoef2('all',c,s,1);
    A1 = appcoef2(c,s,'haar',1);
    V1img = wcodemat(V1,255,'mat',1);
    H1img = wcodemat(H1,255,'mat',1);
    D1img = wcodemat(D1,255,'mat',1);
    A1img = wcodemat(A1,255,'mat',1);

    [H2,V2,D2] = detcoef2('all',c,s,2);
    A2 = appcoef2(c,s,'haar',2);
    V2img = wcodemat(V2,255,'mat',1);
    H2img = wcodemat(H2,255,'mat',1);
    D2img = wcodemat(D2,255,'mat',1);
    A2img = wcodemat(A2,255,'mat',1);

    A1Coeff.V1img = V1img;
    A1Coeff.H1img = H1img;
    A1Coeff.D1img = D1img;
    A1Coeff.A1img = A1img;
    
    A2Coeff.V2img = V2img;
    A2Coeff.H2img = H2img;
    A2Coeff.D2img = D2img;
    A2Coeff.A2img = A2img;

    subplot(2,2,1);
    imagesc(A1img);
    colormap pink(255);
    title('Approximation Coef. of Level 1');

    subplot(2,2,2);
    imagesc(H1img);
    title('Horizontal detail Coef. of Level 1');

    subplot(2,2,3);
    imagesc(V1img);
    title('Vertical detail Coef. of Level 1');

    subplot(2,2,4);
    imagesc(D1img);
    title('Diagonal detail Coef. of Level 1'); 
end
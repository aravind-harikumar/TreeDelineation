function cmap = getColorMap(c1, c2, c3, c4, size, plotOn)
    % Prepare example data rage 1-3
    if(plotOn)
        Z=peaks;
        Z=peaks./max(Z(:));
        Z=(Z+1)*3/2;
    end

    n1=20; n2=20;
    cmap=[linspace(c1(1),c2(1),n1); linspace(c1(2),c2(2),n1); linspace(c1(3),c2(3),n1)];
    cmap(:,end+1:end+n2)=[linspace(c2(1),c3(1),n2);linspace(c2(2),c3(2),n2);linspace(c2(3),c3(3),n2)];
    cmap(:,end+1:end+n2)=[linspace(c3(1),c4(1),n2);linspace(c3(2),c4(2),n2);linspace(c3(3),c4(3),n2)];
    cmap = cmap';
    %colormap(cmap)

    %Plot
    if(plotOn)
        surf(Z);
        colormap(cmap)
        colorbar;
        caxis([1 3]);
    end
end
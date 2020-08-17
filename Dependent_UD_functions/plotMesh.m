function plotMesh(csm, xDiv, yDiv)
    hold on;
    [XXXX,YYYY] = meshgrid(linspace(0, size(csm,2), xDiv),linspace(1, size(csm,1), yDiv));
    plot(XXXX,YYYY,'k')
    [XXXX,YYYY] = meshgrid(linspace(0, size(csm,1), xDiv),linspace(1, size(csm,2), yDiv));
    plot(YYYY,XXXX,'k')
    title('CSM');
    hold off;
end
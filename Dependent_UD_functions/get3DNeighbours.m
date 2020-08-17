function neighbors = get3DNeighbours(inputMatrix,index, neighSize)
    s=size(inputMatrix);
    N=length(s);
    [c1{1:N}]=ndgrid(1:neighSize);
    c2(1:N)={2};
    offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
    neighbors = inputMatrix(index+size(inputMatrix,1)+offsets);
end
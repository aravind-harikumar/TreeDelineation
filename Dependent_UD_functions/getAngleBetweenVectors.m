function angleBwVextors = getAngleBetweenVectors(vec1, vec2)
    if(sum(vec1) == 0);
        vec1 = vec1 + 0.000000001;
    end
    if(sum(vec2) == 0);
        vec2 = vec2 + 0.000000001;
    end
    angleBwVextors =  acosd((vec1*vec2')/(norm(vec1)*norm(vec2)));
end
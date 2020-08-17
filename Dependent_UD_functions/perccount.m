function  perccount(jj,maxjj)
    pc_done  =  round((jj/maxjj)*100, 0);
    if(pc_done>10)
        fprintf(1,'\b\b\b');
    else
        fprintf(1,'\b\b');
    end
    fprintf(1,'%s%%',num2str(pc_done));
end

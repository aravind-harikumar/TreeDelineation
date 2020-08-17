function [csm3D,cdm3D] = get3Dcsmcsm(csm, cdm, zDiv)      
    csm = repmat(csm,[1 1 zDiv]);
    csm3D = csm(:);
    
    cdm = repmat(cdm,[1 1 zDiv]);
    cdm3D = cdm(:); 
end
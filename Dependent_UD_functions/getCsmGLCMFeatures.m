function ImageFeatureArr = getCsmGLCMFeatures(csm, f2D, zDiv, plotOn2d)
        
        f2D.autoc = zeros(size(csm));
        f2D.contr = zeros(size(csm));
        f2D.corrm = zeros(size(csm));
        f2D.corrp = zeros(size(csm));
        f2D.cprom = zeros(size(csm));
        f2D.cshad = zeros(size(csm));
        f2D.dissi = zeros(size(csm));
        f2D.energ = zeros(size(csm));
        f2D.entro = zeros(size(csm));
        f2D.homom = zeros(size(csm));
        f2D.homop = zeros(size(csm));
        f2D.maxpr = zeros(size(csm));
        f2D.sosvh = zeros(size(csm));
        f2D.savgh = zeros(size(csm));
        f2D.svarh = zeros(size(csm));
        f2D.senth = zeros(size(csm));
        f2D.dvarh = zeros(size(csm));
        f2D.denth = zeros(size(csm));
        f2D.inf1h = zeros(size(csm));
        f2D.inf2h = zeros(size(csm));
        f2D.indnc = zeros(size(csm));
        f2D.idmnc = zeros(size(csm));

        neighSize = 3;
        for r = (neighSize-1):1:size(csm,1)-(neighSize-2)
            for c = (neighSize-1):1:size(csm,2)-(neighSize-2)
                
                indx =  sub2ind(size(csm),r,c)-size(csm,1);
                neighbors = getNeighbours(csm,indx,neighSize);
                
                %neighbors = ((neighbors-min(neighbors(:)))/( max(neighbors(:)) - min(neighbors(:)) ))*255;
                neighbors = double(int32(neighbors));
                %rshMat  =reshape(neighbors,[3 3]);
                GLCM2 = graycomatrix(neighbors,'Offset',[0 1; -1 0; -1 -1; -1 1;],'Symmetric', false, 'NumLevels', 8, 'GrayLimits', [0 255]); % -1 0; -1 -1; -1 1
                %glcmfeatures = GLCM_Features1(GLCM2,0);
                glcmfeatures = harFeatures(GLCM2, 12);
                glcmfeatures(isnan(glcmfeatures))=0;
                
                %glcmfeatures = graycoprops(GLCM2,{'contrast','homogeneity'});

        
%                 glcmfeatures.autoc(isnan(glcmfeatures.autoc))=0;
%                 glcmfeatures.contr(isnan(glcmfeatures.contr))=0;
%                 glcmfeatures.corrm(isnan(glcmfeatures.corrm))=0;
%                 glcmfeatures.corrp(isnan(glcmfeatures.corrp))=0;        
%                 glcmfeatures.cprom(isnan(glcmfeatures.cprom))=0;
%                 glcmfeatures.cshad(isnan(glcmfeatures.cshad))=0;
%                 glcmfeatures.dissi(isnan(glcmfeatures.dissi))=0;
%                 glcmfeatures.energ(isnan(glcmfeatures.energ))=0;
%                 glcmfeatures.entro(isnan(glcmfeatures.entro))=0;
%                 glcmfeatures.homom(isnan(glcmfeatures.homom))=0;
%                 glcmfeatures.homop(isnan(glcmfeatures.homop))=0;
%                 glcmfeatures.maxpr(isnan(glcmfeatures.maxpr))=0;
%                 glcmfeatures.sosvh(isnan(glcmfeatures.sosvh))=0;
%                 glcmfeatures.savgh(isnan(glcmfeatures.savgh))=0;
%                 glcmfeatures.svarh(isnan(glcmfeatures.svarh))=0;
%                 glcmfeatures.senth(isnan(glcmfeatures.senth))=0;
%                 glcmfeatures.dvarh(isnan(glcmfeatures.dvarh))=0;
%                 glcmfeatures.denth(isnan(glcmfeatures.denth))=0;
%                 glcmfeatures.inf1h(isnan(glcmfeatures.inf1h))=0;
%                 glcmfeatures.inf2h(isnan(glcmfeatures.inf2h))=0;
%                 glcmfeatures.indnc(isnan(glcmfeatures.indnc))=0;
%                 glcmfeatures.idmnc(isnan(glcmfeatures.idmnc))=0;
                
                f2D.autoc(r,c) = mean(glcmfeatures(:,1));
                f2D.contr(r,c) = mean(glcmfeatures(:,2));
                f2D.corrm(r,c) = mean(glcmfeatures(:,3));
                f2D.corrp(r,c) = mean(glcmfeatures(:,4));
                f2D.cprom(r,c) = mean(glcmfeatures(:,5));
                f2D.cshad(r,c) = mean(glcmfeatures(:,6));
                f2D.dissi(r,c) = mean(glcmfeatures(:,7));
                f2D.energ(r,c) = mean(glcmfeatures(:,8));
                f2D.entro(r,c) = mean(glcmfeatures(:,9));
                f2D.homom(r,c) = mean(glcmfeatures(:,10));
                f2D.homop(r,c) = mean(glcmfeatures(:,11));
                f2D.maxpr(r,c) = mean(glcmfeatures(:,12));
%                 f2D.sosvh(r,c) = mean(glcmfeatures.sosvh);
%                 f2D.savgh(r,c) = mean(glcmfeatures.savgh);
%                 f2D.svarh(r,c) = mean(glcmfeatures.svarh);
%                 f2D.senth(r,c) = mean(glcmfeatures.senth);
%                 f2D.dvarh(r,c) = mean(glcmfeatures.dvarh);
%                 f2D.denth(r,c) = mean(glcmfeatures.denth);
%                 f2D.inf1h(r,c) = mean(glcmfeatures.inf1h);
%                 f2D.inf2h(r,c) = mean(glcmfeatures.inf2h);
%                 f2D.indnc(r,c) = mean(glcmfeatures.indnc);
%                 f2D.idmnc(r,c) = mean(glcmfeatures.idmnc);

                %fImg(r,c) = dd(1);

            end
        end
        
        if(plotOn2d)
        
            f3=figure('name','CSM Texture Features'); % Next Figure
            set(f3, 'Position', [10 25 950 500]);

            subplot(3,4,1);
            factor =0.01;
            f2D.autoc = imgaussfilt(f2D.autoc,factor);
            f2D.autoc = normalize(f2D.autoc);
            imagesc(f2D.autoc);%1
            title('autoc');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,2);
            f2D.contr = imgaussfilt(f2D.contr,factor);
            f2D.contr = normalize(f2D.contr);
            imagesc(f2D.contr);%2
            title('contr');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,3);
            f2D.corrm = imgaussfilt(f2D.corrm,factor);
            f2D.corrm = normalize(f2D.corrm);
            imagesc(f2D.corrm);%3
            title('corrm');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,4);
            f2D.corrp = imgaussfilt(f2D.corrp,factor);
            f2D.corrp = normalize(f2D.corrp);
            imagesc(f2D.corrp);%4
            title('corrp');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,5);
            f2D.cprom = imgaussfilt(f2D.cprom,factor);
            f2D.cprom = normalize(f2D.cprom);
            imagesc(f2D.cprom);%5
            title('cprom');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,6);   
            f2D.cshad = imgaussfilt(f2D.cshad,factor);
            f2D.cshad = normalize(f2D.cshad);
            imagesc(f2D.cshad);%6
            title('cshad');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,7);
            f2D.dissi = imgaussfilt(f2D.dissi,factor);
            f2D.dissi = normalize(f2D.dissi);
            imagesc(f2D.dissi);%7
            title('dissi');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,8);
            f2D.energ = imgaussfilt(f2D.energ,factor);
            f2D.energ = normalize(f2D.energ);
            imagesc(f2D.energ);%8
            title('energ');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,9);
            f2D.entro = imgaussfilt(f2D.entro,factor);
            f2D.entro = normalize(f2D.entro);
            imagesc(f2D.entro);
            title('entro');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,10);
            f2D.homom = imgaussfilt(f2D.homom,factor);
            f2D.homom = normalize(f2D.homom);
            imagesc(f2D.homom);
            title('homom');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,11);
            f2D.homop = imgaussfilt(f2D.homop,factor);
            f2D.homop = normalize(f2D.homop);
            imagesc(f2D.homop);
            title('homop');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            subplot(3,4,12);
            f2D.maxpr = imgaussfilt(f2D.maxpr,factor);
            f2D.maxpr = normalize(f2D.maxpr);
            imagesc(f2D.maxpr);
            title('maxpr');
            set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

            [~,h3]=suplabel('2D TEXTURE FEATURES' ,'t');
            set(h3,'FontSize',18);
        end
        % ---------------------------------------------------------------%  
        % 2D to 3D conversion to allocate to voxels (repeat layers N times)       
        % ---------------------------------------------------------------%  
        f3D.filtcsm = repmat(f2D.filtcsmArr{1},[1 1 zDiv]);
        f3D.filtcsm = f3D.filtcsm(:);
        
        f3D.entrImage = repmat(f2D.entrImage,[1 1 zDiv]);
        f3D.entrImage = f3D.entrImage(:);
        
        f3D.J = repmat(f2D.J,[1 1 zDiv]);
        f3D.J = f3D.J(:);
        
        f3D.diff1 = repmat(f2D.diffFetArr{1},[1 1 zDiv]);
        f3D.diff1 = f3D.diff1(:);
        
        f3D.diff2 = repmat(f2D.diffFetArr{2},[1 1 zDiv]);
        f3D.diff2 = f3D.diff2(:);
        
        f3D.diff3 = repmat(f2D.diffFetArr{3},[1 1 zDiv]);
        f3D.diff3 = f3D.diff3(:);
        
        f3D.autoc = repmat(f2D.autoc,[1 1 zDiv]);
        f3D.autoc = f3D.autoc(:);
        
        f3D.contr = repmat(f2D.contr,[1 1 zDiv]);
        f3D.contr = f3D.contr(:);
        
        f3D.corrm = repmat(f2D.corrm,[1 1 zDiv]);
        f3D.corrm = f3D.corrm(:);
        
        f3D.corrp = repmat(f2D.corrp,[1 1 zDiv]);
        f3D.corrp = f3D.corrp(:);
        
        f3D.cprom = repmat(f2D.cprom,[1 1 zDiv]);
        f3D.cprom = f3D.cprom(:);
        
        f3D.cshad = repmat(f2D.cshad,[1 1 zDiv]);
        f3D.cshad = f3D.cshad(:);
        
        f3D.dissi = repmat(f2D.dissi,[1 1 zDiv]);
        f3D.dissi = f3D.dissi(:);
        
        f3D.energ = repmat(f2D.energ,[1 1 zDiv]);
        f3D.energ = f3D.energ(:);
    
        f3D.entro = repmat(f2D.entro,[1 1 zDiv]);
        f3D.entro = f3D.entro(:);
        
        f3D.homom = repmat(f2D.homom,[1 1 zDiv]);
        f3D.homom = f3D.homom(:);
                
        f3D.homop = repmat(f2D.homop,[1 1 zDiv]);
        f3D.homop = f3D.homop(:);
        
        f3D.maxpr = repmat(f2D.maxpr,[1 1 zDiv]);
        f3D.maxpr = f3D.maxpr(:);
        
        ImageFeatureArr = [f3D.filtcsm f3D.diff1 f3D.diff2 f3D.entrImage f3D.J]; %5
        ImageFeatureArr = [ImageFeatureArr f3D.autoc f3D.contr f3D.corrm f3D.corrp f3D.cprom f3D.cshad]; %6
        ImageFeatureArr = [ImageFeatureArr f3D.dissi f3D.energ f3D.entro f3D.homom f3D.homop f3D.maxpr]; %7
        
        
        ImageFeatureArr = [ImageFeatureArr f3D.autoc f3D.contr f3D.corrm f3D.corrp f3D.cprom f3D.cshad f3D.dissi...
        f3D.energ f3D.entro f3D.homom f3D.homop f3D.maxpr]; 
                
        
        %ImageFeatureArr = [ImageFeatureArr f3D.sosvh f3D.savgh f3D.svarh f3D.senth f3D.dvarh f3D.denth f3D.inf1h]; %6
        %ImageFeatureArr = [ImageFeatureArr f3D.inf2h f3D.homom f3D.indnc f3D.idmnc]; %4

        %[ax1,h1]=suplabel('SLICES (distance to stem)'); [ax2,h2]=suplabel('3D GLCM TEXTURE FEATURES','y'); [ax3,h3]=suplabel('3D TEXTURE FEATURES' ,'t');
        %set(h1,'FontSize',18); set(h2,'FontSize',18); set(h3,'FontSize',20)
end


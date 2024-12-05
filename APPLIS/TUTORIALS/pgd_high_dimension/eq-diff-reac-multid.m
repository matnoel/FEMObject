% Applications numeriques pour la publication "Tensor-based methods and
% Proper Generalized Decompositions for the numerical solution of high
% dimensional problems: alternative definitions", G. Bonithon, A. Nouy,
% 2011.

% Exemple 1 :
% Equation de diffusion-reaction -\Delta u + \alpha(x_1,...,x_n)u = f
% Analyse du comportement de la G-PGD et de l'impact des updates par rapport
% a l'elevation en dimension, l'ellipticite, et la tangentialite de
% l'operateur

%% Preambule

clear();
bPrint = false;
bElliptic = false;
bTangential = false;
bvNSTUpSolsComputed = zeros(1,4);
bvNUpSolsComputed = zeros(1,1000);
bManufactRef = true;
strErrCompMode = 'reference';
dFreq = 10;
dAmp = 30;
iDim = 3;
iNbElem = 30;
if (bElliptic)
    strEllipMode = 'Elliptic';
else
    strEllipMode = 'Non elliptic';
end
if (bTangential)
    strTangentialMode = 'tangential';
else
    strTangentialMode = 'non tangential';
end
strDim = [num2str(iDim),'D-'];
strPath = [pwd,'/Publi/Ex1/'];

% Generation du modele
domD1d = DOMAIN(1);
modM1d = mesh(domD1d,iNbElem);
modM1d = createddlnode(modM1d,DDL('u'));
modM1d = addcl(modM1d,POINT([0;1]),'u',0);

% Terme de diffusion (Id\otimes\Delta) et terme source pour une dimension
% Formulation continue
blinfD00 = BILINFORM(0,0);
blinfD11 = BILINFORM(1,1);
dvCoord = getcoord(getnode(modM1d));
dvSrc = ones(size(dvCoord,1),1);
linfSrc = LINFORM(0,dvSrc,0);

% Formulation discrete
dmD00 = blinfD00{modM1d}(:,:);
dmD11 = blinfD11{modM1d}(:,:);
dvSrc1d = linfSrc{modM1d}{:};

% Assemblage des dimensions
cvdvSrc(1:iDim) = {dvSrc1d};
cmdmSys(1:iDim,1:iDim) = {dmD00};
cmdmSys(logical(eye(iDim))) = {dmD11};

% Terme de reaction
% Formulation continue
dvWeight = dAmp*cos(dFreq*pi*dvCoord(:,1));
% dvWeight = ones(length(dvCoord),1);
if (bElliptic)
    dvWeight = abs(dvWeight);
end
blinfD00 = BILINFORM(0,0,dvWeight,0);
if (bTangential)
    blinfD00one = BILINFORM(0,0);
end

% Formulation discrete
dmD00 = blinfD00{modM1d}(:,:);
if (bTangential)
    dmD00one = blinfD00one{modM1d}(:,:);
end

% Assemblage des dimensions
if (bTangential)
    cmdmSys(iDim+1:2*iDim,1:iDim) = {dmD00one};
    cmdmSys(logical([zeros(iDim);eye(iDim)])) = {dmD00};
else
    cmdmSys(iDim+1,1:iDim) = {dmD00};
end

% Passage aux sepmatrix
sepmSys = SEPMATRIX(cmdmSys);
sepvSrc = SEPMATRIX(cvdvSrc);

% Solveur par defaut
solvDefault = SEPSOLVER(getdim(sepmSys),'adjoint',0,'inittype','one',...
    'residual',0,'maxorder',20,'tol',1e-6,'update',0,....
    'errorindicator',strErrCompMode,'ortho',0,'updateeps',0,...
    'itercrit',5e-3);
if (bElliptic || bTangential)
    solvDefault = setparam(solvDefault,'maxorder',100);
end
cvstrMarkers = {'none','+','*','o','^','s','p','x',...
    'none','+','*','o','^','s','p','x'};

solvCurrent = setparam(solvDefault,'maxorder',12,'itercrit',5e-8,...
    'maxiter',100,'tol',1e-7,'iErrCompStep',1,'inittype','rand',...
    'itercritupdate',1e-10);


%% Influence du nombre d'updates pour le O-Update

solvCurrent = setparam(solvCurrent,'maxorder',10);
figure(iDim);
% iNbUp = min(iDim+1,length(cvstrMarkers)-1);
iNbUp = 0;
bvNUpSolsComputed(1:iNbUp+1)=0;

if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,iNbUp+1);
    cvsepvNUpSols = cell(1,iNbUp);
    cvstNUpResults = cell(1,iNbUp);
end
for i=0:iNbUp
    if (~bvNUpSolsComputed(i+1))
        solvNUp = setparam(solvCurrent,'update',i);
        [cvsepvNUpSols{i+1},cvstNUpResults{i+1}] = solve(sepmSys,sepvSrc,solvNUp);
        bvNUpSolsComputed(i+1) = true;
    end
    if (~i)
        semilogy(cvstNUpResults{i+1}.error,'r','Marker',cvstrMarkers{i+1});
    else
        semilogy(cvstNUpResults{i+1}.error,'Marker',cvstrMarkers{i+1});
    end
    hold on;
    if (~i)
        cvstrLegends{i+1} = 'G-PGD';
        if (~bvNSTUpSolsComputed(1))
            cvsepvNSTUpSols{1} = cvsepvNUpSols{i+1};
            cvstNSTUpResults{1} = cvstNUpResults{i+1};
            bvNSTUpSolsComputed(1) = true;
        end
    elseif (i == iNbUp)
        cvstrLegends{i+1} = ['(',num2str(i),')G-PGD'];
        if (~bvNSTUpSolsComputed(2))
            cvsepvNSTUpSols{2} = cvsepvNUpSols{i+1};
            cvstNSTUpResults{2} = cvstNUpResults{i+1};
            bvNSTUpSolsComputed(2) = true;
        end
    else
        cvstrLegends{i+1} = ['(',num2str(i),')G-PGD'];
    end
end

xlabel('Rank');
ylabel('Error');
legend(cvstrLegends);

strTitle = [strDim,'GPGD convergence (',strErrCompMode,' error) - ',...
    strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

if bPrint
    bPrint = false;
    title('');
    strFileName = lower(strrep([strDim,'nup-GPGD-cv-',strErrCompMode,'-err-',...
        strEllipMode,'-',strTangentialMode],' ','-'));
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end

%% Calcul d'une solution de reference : solution "manufacturee" ou FEM en 2D et 3D et MSVD(GPGD) au dela

if (bManufactRef)
    % Construction d'une solution sur la base des polynomes : accessible pour la PGD
    % sans etre trivial comme les fonctions propres du laplacien
    iNbModes = 20;
    dvOpenCoord = dvCoord(2:end-1); % suppression des noeuds au bord, mieux?
    iDeg = 10;
    imPowers = randi(iDeg+1,iNbModes,iDim)-1;
    dvAlpha = ones(1,iNbModes);
    cmdvModes = cell(iNbModes,iDim);
    for i=1:iNbModes
        for j=1:iDim
            cmdvModes{i,j} = dvOpenCoord.*(1-dvOpenCoord).*dvOpenCoord.^imPowers(i,j);
            dvAlpha(i) = dvAlpha(i)*norm(cmdvModes{i,j});
            cmdvModes{i,j} = cmdvModes{i,j}/norm(cmdvModes{i,j});
        end
    end
    sepvRefSol = SEPMATRIX(cmdvModes,dvAlpha);
    sepvSrc = sepmSys*sepvRefSol;
    % sepvSrc = sepvRefSol;
    
    % Solution purement aleatoire, le plus simple pour avoir des
    % "vraies" courbes de convergence
    %     iNbModes = 20;
    %     cmdvModes = cell(iNbModes,iDim);
    %     dvOpenCoord = dvCoord(2:end-1); % suppression des noeuds au bord, mieux?
    %     iNbIntPts = length(dvOpenCoord);
    %     dvAlpha = ones(1,iNbModes);
    %     iNbRandPts = 20;
    %     for i=1:iNbModes
    %         for j=1:iDim
    %             %dvRandElemTens = rand(iNbIntPts,1);
    %             dvRandElemTens = interp1(0:1/(iNbRandPts-1):1,rand(1,iNbRandPts),dvOpenCoord,'spline');
    %             cmdvModes{i,j} = dvOpenCoord.*(1-dvOpenCoord).*dvRandElemTens;
    %             dvAlpha(i) = dvAlpha(i)*norm(cmdvModes{i,j});
    %             cmdvModes{i,j} = cmdvModes{i,j}/norm(cmdvModes{i,j});
    %         end
    %     end
    %     plot(dvOpenCoord,cmdvModes{2,3})
    %     sepvRefSol = SEPMATRIX(cmdvModes,dvAlpha);
    %     sepvSrc = sepmSys*sepvRefSol;
    
    % Construction d'une solution sur la base des vecteurs propres du
    % laplacien : trop "fleche", conduit a des courbes de convergence
    % triviales (PGD=U-PGD=SVD)
    %     iNbModes = 10;
    %     dvRefSeq = 1./((1:iNbModes).^3);
    %     %dvAlpha = rand(1,iNbModes).*dvRefSeq;
    %     dvAlpha = dvRefSeq;
    %     %dvAlpha = ones(1,iNbModes);
    %     iMaxFreq = 30;
    %     dmFreq = pi*randi(iMaxFreq,iNbModes,iDim);
    %     cmdvModes = cell(iNbModes,iDim);
    %     dvOpenCoord = dvCoord(2:end-1); % suppression des noeuds au bord, mieux?
    %     for i=1:iNbModes
    %         for j=1:iDim
    %             cmdvModes{i,j} = sin(dmFreq(i,j)*dvOpenCoord)/...
    %                 sqrt(1/2-sin(2*dmFreq(i,j))/(4*dmFreq(i,j)));
    %         end
    %     end
    %     sepvRefSol = SEPMATRIX(cmdvModes,dvAlpha);
    %     sepvSrc = sepmSys*sepvRefSol;
    
    % MSVD du terme source
    %     solvSvd = setparam(solvCurrent,'maxorder',50,'update',3,...
    %         'errorindicator','residual','tol',5e-3);
    %     [sepvSrc,stResult] = multisvd(sepvSrc,solvSvd);
    %     semilogy(stResult.error);
    %     xlabel('Rank');
    %     ylabel('Error');
    %     strTitle = [strDim,'MSVD convergence for the source term'];
    %     title(strTitle);
    
    if (strcmp(strErrCompMode,'reference'))
        solvCurrent = setparam(solvCurrent,'reference',sepvRefSol);
        if (bElliptic)
            solvCurrent = setparam(solvCurrent,'metric',sepmSys);
        end
    end
else
    if (iDim<=3)
        % Generation du modele
        domD = DOMAIN(iDim);
        modM = mesh(domD,iNbElem*ones(1,iDim));
        modM = createddlnode(modM,DDL('u'));
        modM = addcl(modM,[],'u',0);
        dmCoord = getcoord(getnode(modM));
        
        % Formulation continue
        % Terme de diffusion et terme source
        blinfD11 = BILINFORM(1,1);
        dvSrc = ones(size(dmCoord,1),1);
        linfSrc = LINFORM(0,dvSrc,0);
        
        % Terme de reaction
        if (bTangential)
            if (bElliptic)
                dvWeight = sum(abs(dAmp*cos(dFreq*pi*dmCoord)),2);
            else
                dvWeight = sum(dAmp*cos(dFreq*pi*dmCoord),2);
            end
        else
            dvWeight = prod(dAmp*cos(dFreq*pi*dmCoord),2);
            if (bElliptic)
                dvWeight = abs(dvWeight);
            end
        end
        blinfD00 = BILINFORM(0,0,dvWeight,0);
        
        % Formulation discrete
        dmD11 = blinfD11{modM}(:,:);
        dvSrc = linfSrc{modM}(:);
        dmD00 = blinfD00{modM}(:,:);
        cmdmSys = dmD11 + dmD00;
        
        % Resolution et reference pour le calcul d'erreur PGD
        dvSol = solve(cmdmSys,dvSrc);
        dmRefSol = reshape(full(dvSol),(iNbElem-1)*ones(1,iDim));
        
        % Trace et sauvegarde des figures
        figure;
        plot(dvSol,modM);
        strTitle = [strDim,'FEM Solution - ',strEllipMode,' ',strTangentialMode,' problem'];
        title(strTitle);
        if bPrint
            bPrint = false;
            strFileName = lower(strrep([strDim,'FEM-Sol-',strEllipMode,'-',strTangentialMode],' ','-'));
            myprint(strPath,strFileName,'epsc2');
            save([strPath,strFileName],'dmRefSol','dvSol');
        end
        if (strcmp(strErrCompMode,'reference'))
            solvCurrent = setparam(solvCurrent,'reference',dmRefSol);
        end
    else
        iNbNUp = 3;
        solvCurrent = setparam(solvDefault,'errorindicator','residual','update',iNbNUp,...
            'maxorder',100,'itercrit',5e-3,'maxiter',50,'tol',1e-5,'iErrCompStep',1);
        
        [sepvRefSol,result] = solve(sepmSys,sepvSrc,solvCurrent);
        figure;
        semilogy(result.error);
        hold on;
        [sepvRefSol,result] = multisvd(sepvRefSol,solvCurrent);
        semilogy(result.error,'Marker','+');
        
        xlabel('Rank');
        ylabel('Error');
        legend(['(',num2str(iNbNUp),')O-GPGD'],...
            ['((',num2str(iNbNUp),')O-MSVD (',num2str(iNbNUp),')O-GPGD)']);
        strTitle = ['Convergence of a',strDim,'separated refernce solution '...
            '(residual error) - ',strEllipMode,' ',strTangentialMode,' problem'];
        title(strTitle);
        if bPrint
            bPrint = false;
            strFileName = lower(strrep([strDim,'GPGD-MSVD-ref-cv-',strErrCompMode,'-err-',...
                strEllipMode,'-',strTangentialMode],' ','-'));
            save([strPath,strFileName],'sepvRefSol');
            myprint(strPath,strFileName,'epsc2');
        end
        if (strcmp(strErrCompMode,'reference'))
            solvCurrent = setparam(solvCurrent,'reference',sepvRefSol);
        end
    end
end


%% Comparaison des modes d'update

bRes = false;
if (bRes)
    strKindOfPGD = 'MR';
else
    strKindOfPGD = 'G';
end
% bvNSTUpSolsComputed(:)=0;
solvCurrent = setparam(solvCurrent,'maxorder',20,'tol',1e-12,'residual',bRes,...
    'maxiter',100);
if (~exist('cvsepvNSTUpSols','var'))
    cvsepvNSTUpSols = cell(1,4);
    cvstNSTUpResults = cell(1,4);
end

if (~bvNSTUpSolsComputed(1))
    [cvsepvNSTUpSols{1},cvstNSTUpResults{1}] = solve(sepmSys,sepvSrc,solvCurrent);
    bvNSTUpSolsComputed(1) = true;
end
if (~bvNUpSolsComputed(1))
    cvsepvNUpSols{1} = cvsepvNSTUpSols{1};
    cvstNUpResults{1} = cvstNSTUpResults{1};
    bvNUpSolsComputed(1) = true;
end

% bvNSTUpSolsComputed(1:4)=0;
figure;
cvstrLegends{1} = [strKindOfPGD,'-PGD'];
semilogy(cvstNSTUpResults{1}.error,'r');
hold on

if (~bvNSTUpSolsComputed(3))
    solvAUp = setparam(solvCurrent,'alphaupdate',true);
    [cvsepvNSTUpSols{3},cvstNSTUpResults{3}] = solve(sepmSys,sepvSrc,solvAUp);
    bvNSTUpSolsComputed(3) = true;
end

cvstrLegends{2} = ['(S)',strKindOfPGD,'-PGD'];
semilogy(cvstNSTUpResults{3}.error,'Marker',cvstrMarkers{3});

if (~bvNSTUpSolsComputed(2))
    iNbUp = min(iDim+1,length(cvstrMarkers)-1);
    iNbUp = 1;
    solvNUp = setparam(solvCurrent,'update',iNbUp);
    [cvsepvNSTUpSols{2},cvstNSTUpResults{2}] = solve(sepmSys,sepvSrc,solvNUp);
    bvNSTUpSolsComputed(2) = true;
end
if (~bvNUpSolsComputed(iNbUp+1))
    cvsepvNUpSols{iNbUp+1} = cvsepvNSTUpSols{2};
    cvstNUpResults{iNbUp+1} = cvstNSTUpResults{2};
    bvNUpSolsComputed(iNbUp+1) = true;
end
cvstrLegends{3} = ['(',num2str(iNbUp),')',strKindOfPGD,'-PGD'];
semilogy(cvstNSTUpResults{2}.error,'Marker',cvstrMarkers{2});

if (iDim<3)
    if (~bvNSTUpSolsComputed(4))
        solvAUp = setparam(solvCurrent,'updatetucker',true);
        [cvsepvNSTUpSols{4},cvstNSTUpResults{4}] = solve(sepmSys,sepvSrc,solvAUp);
        bvNSTUpSolsComputed(4) = true;
    end
    cvstrLegends{4} = ['(T)',strKindOfPGD,'-PGD'];
    semilogy(cvstNSTUpResults{4}.error,'Marker',cvstrMarkers{4});
end

xlabel('Rank');
ylabel('Error');
ylim([0,100]);
if (iDim<3)
    cvstrLegends = cvstrLegends(1:4);
else
    cvstrLegends = cvstrLegends(1:3);
end
legend(cvstrLegends,'Location','NorthWest');
strTitle = [strDim,strKindOfPGD,'-PGD convergence (',strErrCompMode,' error) - ',...
    strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);
if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,'nstup-',strKindOfPGD,'PGD-cv-',strErrCompMode,'-err-',...
        strEllipMode,'-',strTangentialMode],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end


%% Erreur vs dimension pour un rang fixe

iMaxOrder = 10;
iMaxDim = 50;
iNbSols = iMaxDim-1;
if (~exist('dvFixedRankErr0','var'))
    iOldMaxDim = 1;
    cvstrLegends = cell(1,2);
    cvstrLegends{1} = 'G-PGD';
    cvstrLegends{2} = '(1)G-PGD';
    cvstrLegends{3} = '(2)G-PGD';
    sepmSys = cell(1,iNbSols);
    sepvSrc = cell(1,iNbSols);
    sepvRefSol = cell(1,iNbSols);
    sepvFixedRankSol0 = cell(1,iNbSols);
    sepvFixedRankSol1 = cell(1,iNbSols);
    sepvFixedRankSol2 = cell(1,iNbSols);
    stFixedRankRes0 = cell(1,iNbSols);
    stFixedRankRes1 = cell(1,iNbSols);
    stFixedRankRes2 = cell(1,iNbSols);
    dvFixedRankErr0 = ones(1,iNbSols);
    dvFixedRankErr1 = ones(1,iNbSols);
    dvFixedRankErr2 = ones(1,iNbSols);
else
    iOldNbSols = length(dvFixedRankErr0);
    iOldMaxDim = iOldNbSols+1;
    iNbNewDims = iNbSols-iOldNbSols;
    sepmSys = [sepmSys,cell(1,iNbNewDims)];
    sepvSrc = [sepvSrc,cell(1,iNbNewDims)];
    sepvRefSol = [sepvRefSol,cell(1,iNbNewDims)];
    sepvFixedRankSol0 = [sepvFixedRankSol0,cell(1,iNbNewDims)];
    sepvFixedRankSol1 = [sepvFixedRankSol1,cell(1,iNbNewDims)];
    sepvFixedRankSol2 = [sepvFixedRankSol2,cell(1,iNbNewDims)];
    stFixedRankRes0 = [stFixedRankRes0,cell(1,iNbNewDims)];
    stFixedRankRes1 = [stFixedRankRes1,cell(1,iNbNewDims)];
    stFixedRankRes2 = [stFixedRankRes2,cell(1,iNbNewDims)];
    dvFixedRankErr0 = [dvFixedRankErr0,ones(1,iNbNewDims)];
    dvFixedRankErr1 = [dvFixedRankErr1,ones(1,iNbNewDims)];
    dvFixedRankErr2 = [dvFixedRankErr2,ones(1,iNbNewDims)];
end

iNbElem = 1000;
iCompDim = 0;
% ivCompDims = zeros(1,1);
for iDim=iOldMaxDim+1:iMaxDim
    if (iDim<10 || ~rem(iDim,10))
        iCompDim = iCompDim + 1;
        ivCompDims(iCompDim) = iDim;
        % Generation du modele
        domD1d = DOMAIN(1);
        modM1d = mesh(domD1d,iNbElem);
        modM1d = createddlnode(modM1d,DDL('u'));
        modM1d = addcl(modM1d,POINT([0;1]),'u',0);
        
        % Terme de diffusion (Id\otimes\Delta) et terme source pour une dimension
        % Formulation continue
        blinfD00 = BILINFORM(0,0);
        blinfD11 = BILINFORM(1,1);
        dvCoord = getcoord(getnode(modM1d));
        dvSrc = ones(size(dvCoord,1),1);
        linfSrc = LINFORM(0,dvSrc,0);
        
        % Formulation discrete
        dmD00 = blinfD00{modM1d}(:,:);
        dmD11 = blinfD11{modM1d}(:,:);
        dvSrc1d = linfSrc{modM1d}{:};
        
        % Assemblage des dimensions
        % iDim = 10;
        cvdvSrc = cell(1,iDim);
        cmdmSys = cell(iDim,iDim);
        for i=1:iDim
            cvdvSrc{i} = dvSrc1d;
            for j=1:iDim
                if (i==j)
                    cmdmSys{i,j} = dmD11;
                else
                    cmdmSys{i,j} = dmD00;
                end
            end
        end
        
        % Terme de reaction
        % Formulation continue
        dvWeight = dAmp*cos(dFreq*pi*dvCoord(:,1));
        % dvWeight = dAmp*exp(-dFreq*pi*dvCoord(:,1).^2);
        % dvWeight = ones(length(dvCoord),1);
        if (bElliptic)
            dvWeight = abs(dvWeight);
        end
        blinfD00 = BILINFORM(0,0,dvWeight,0);
        if (bTangential)
            blinfD00one = BILINFORM(0,0);
        end
        
        % Formulation discrete
        dmD00 = blinfD00{modM1d}(:,:);
        if (bTangential)
            dmD00one = blinfD00one{modM1d}(:,:);
        end
        
        % Assemblage des dimensions
        if (bTangential)
            cmdmSys = [cmdmSys;cell(iDim,iDim)];
            cvdmSys = cell(1,iDim);
            for i=1:iDim-1
                cvdmSys{i} = dmD00one;
            end
            cvdmSys{end} = dmD00;
            for i=1:iDim
                cmdmSys(iDim+i,:) = circshift(cvdmSys,[1,i]);
            end
        else
            cmdmSys = [cmdmSys;cell(1,iDim)];
            for i=1:iDim
                cmdmSys{iDim+1,i} = dmD00;
            end
        end
        
        % Passage aux sepmatrix
        iSolNum = iDim-1;
        sepmSys{iSolNum} = SEPMATRIX(cmdmSys);
        cvsolvFixedRank0 = SEPSOLVER(getdim(sepmSys{iSolNum}),'adjoint',0,...
            'residual',0,'errorindicator',strErrCompMode,'ortho',0,'updateeps',0,...
            'itercrit',5e-8,'maxiter',100,'tol',1e-7,'iErrCompStep',1,'inittype','rand',...
            'itercritupdate',1e-10,'maxorder',iMaxOrder);
        cvsolvFixedRank1 = setparam(cvsolvFixedRank0,'update',1);
        cvsolvFixedRank2 = setparam(cvsolvFixedRank0,'update',2);
        
        if (bManufactRef)
            % Construction d'une solution sur la base des polynomes : accessible pour la PGD
            % sans etre trivial comme les fonctions propres du laplacien
            iNbModes = 20;
            dvOpenCoord = dvCoord(2:end-1); % suppression des noeuds au bord, mieux?
            iDeg = 10;
            imPowers = randi(iDeg+1,iNbModes,iDim)-1;
            dvAlpha = ones(1,iNbModes);
            cmdvModes = cell(iNbModes,iDim);
            for i=1:iNbModes
                for j=1:iDim
                    cmdvModes{i,j} = dvOpenCoord.*(1-dvOpenCoord).*dvOpenCoord.^imPowers(i,j);
                    dvAlpha(i) = dvAlpha(i)*norm(cmdvModes{i,j});
                    cmdvModes{i,j} = cmdvModes{i,j}/norm(cmdvModes{i,j});
                end
            end
            sepvRefSol{iSolNum} = SEPMATRIX(cmdvModes,dvAlpha);
            sepvSrc{iSolNum} = sepmSys{iSolNum}*sepvRefSol{iSolNum};
            if (strcmp(strErrCompMode,'reference'))
                cvsolvFixedRank0 = setparam(cvsolvFixedRank0,'reference',...
                    sepvRefSol{iSolNum},'metric',sepmSys{iSolNum});
                cvsolvFixedRank1 = setparam(cvsolvFixedRank1,'reference',...
                    sepvRefSol{iSolNum},'metric',sepmSys{iSolNum});
                cvsolvFixedRank2 = setparam(cvsolvFixedRank2,'reference',...
                    sepvRefSol{iSolNum},'metric',sepmSys{iSolNum});
            end
        else
            sepvSrc{iSolNum} = SEPMATRIX(cvdvSrc);
        end
        
        [sepvFixedRankSol0{iSolNum},stFixedRankRes0{iSolNum}] =...
            solve(sepmSys{iSolNum},sepvSrc{iSolNum},cvsolvFixedRank0);
        dvFixedRankErr0(iSolNum) = stFixedRankRes0{iSolNum}.error(iMaxOrder);
        [sepvFixedRankSol1{iSolNum},stFixedRankRes1{iSolNum}] =...
            solve(sepmSys{iSolNum},sepvSrc{iSolNum},cvsolvFixedRank1);
        dvFixedRankErr1(iSolNum) = stFixedRankRes1{iSolNum}.error(iMaxOrder);
        [sepvFixedRankSol2{iSolNum},stFixedRankRes2{iSolNum}] =...
            solve(sepmSys{iSolNum},sepvSrc{iSolNum},cvsolvFixedRank2);
        dvFixedRankErr2(iSolNum) = stFixedRankRes2{iSolNum}.error(iMaxOrder);
    end
end

figure(iMaxOrder);
semilogy(ivCompDims,dvFixedRankErr0(ivCompDims),'r');
hold on
semilogy(ivCompDims,dvFixedRankErr1(ivCompDims),'Marker','+');
hold on
semilogy(ivCompDims,dvFixedRankErr2(ivCompDims),'Marker','*');
xlabel('Dimension');
ylabel('Error');
legend(cvstrLegends,'Location','SouthEast');
strTitle = ['GPGD convergence at rank ',num2str(iMaxOrder),' vs dimension (',...
    strErrCompMode,' error) - ',strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

if bPrint
    bPrint = false;
    strFileName = lower(strrep(['gpgd-cv-rank',num2str(iMaxOrder),'-dim2-',num2str(iMaxDim),'-',...
        strErrCompMode,'-err-',strEllipMode,'-',strTangentialMode],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end


%% Rang vs dimension pour une erreur fixe

iUnits = 1;
iExp = 3;
dMaxErr = iUnits*10^(-iExp);
iMaxDim = 5;
iNbUp = 2;
iNbSols = iMaxDim-1;
if (~exist('cvdvFixedErrRank','var'))
    iOldMaxDim = 1;
    cvsepmSys = cell(1,iNbSols);
    cvsepvSrc = cell(1,iNbSols);
    cvsepvRefSol = cell(1,iNbSols);
    cvcvsepvFixedErrSol(1:iNbUp) = {cell(1,iNbSols)};
    cvcvstFixedErrRes(1:iNbUp) = {cell(1,iNbSols)};
    cvdvFixedErrRank(1:iNbUp) = {ones(1,iNbSols)};
    cvstrLegends = cell(1,iNbUp+1);
    cvstrLegends{1} = 'G-PGD';
    for i=1:iNbUp
        cvstrLegends{i+1} = ['(',num2str(i),')G-PGD'];
    end
else
    iOldNbSols = length(cvdvFixedErrRank{1});
    iOldMaxDim = iOldNbSols+1;
    iNbNewDims = iNbSols-iOldNbSols;
    cvsepmSys = [cvsepmSys,cell(1,iNbNewDims)];
    cvsepvSrc = [cvsepvSrc,cell(1,iNbNewDims)];
    cvsepvRefSol = [cvsepvRefSol,cell(1,iNbNewDims)];
    for i=1:iNbUp+1
        cvcvsepvFixedErrSol{i} = [cvcvsepvFixedErrSol{i},cell(1,iNbNewDims)];
        cvcvstFixedErrRes{i} = [cvcvstFixedErrRes{i},cell(1,iNbNewDims)];
        cvdvFixedErrRank{i} = [cvdvFixedErrRank{i},ones(1,iNbNewDims)];
    end
end

for iDim=iOldMaxDim+1:iMaxDim
    % Terme de diffusion (Id\otimes\Delta) et terme source pour une dimension
    % Formulation continue
    blinfD00 = BILINFORM(0,0);
    blinfD11 = BILINFORM(1,1);
    dvCoord = getcoord(getnode(modM1d));
    dvSrc = ones(size(dvCoord,1),1);
    linfSrc = LINFORM(0,dvSrc,0);
    
    % Formulation discrete
    dmD00 = blinfD00{modM1d}(:,:);
    dmD11 = blinfD11{modM1d}(:,:);
    dvSrc1d = linfSrc{modM1d}{:};
    
    % Assemblage des dimensions
    % iDim = 10;
    cvdvSrc = cell(1,iDim);
    cmdmSys = cell(iDim,iDim);
    for i=1:iDim
        cvdvSrc{i} = dvSrc1d;
        for j=1:iDim
            if (i==j)
                cmdmSys{i,j} = dmD11;
            else
                cmdmSys{i,j} = dmD00;
            end
        end
    end
    
    % Terme de reaction
    % Formulation continue
    dvWeight = dAmp*cos(dFreq*pi*dvCoord(:,1));
    % dvWeight = dAmp*exp(-dFreq*pi*dvCoord(:,1).^2);
    % dvWeight = ones(length(dvCoord),1);
    if (bElliptic)
        dvWeight = abs(dvWeight);
    end
    blinfD00 = BILINFORM(0,0,dvWeight,0);
    if (bTangential)
        blinfD00one = BILINFORM(0,0);
    end
    
    % Formulation discrete
    dmD00 = blinfD00{modM1d}(:,:);
    if (bTangential)
        dmD00one = blinfD00one{modM1d}(:,:);
    end
    
    % Assemblage des dimensions
    if (bTangential)
        cmdmSys = [cmdmSys;cell(iDim,iDim)];
        cvdmSys = cell(1,iDim);
        for i=1:iDim-1
            cvdmSys{i} = dmD00one;
        end
        cvdmSys{end} = dmD00;
        for i=1:iDim
            cmdmSys(iDim+i,:) = circshift(cvdmSys,[1,i]);
        end
    else
        cmdmSys = [cmdmSys;cell(1,iDim)];
        for i=1:iDim
            cmdmSys{iDim+1,i} = dmD00;
        end
    end
    
    % Passage aux sepmatrix
    iSolNum = iDim-1;
    cvsepmSys{iSolNum} = SEPMATRIX(cmdmSys);
    cvsolvFixedRank{1} = SEPSOLVER(getdim(cvsepmSys{iSolNum}),'adjoint',0,...
        'residual',0,'errorindicator',strErrCompMode,'ortho',0,'updateeps',0,...
        'itercrit',5e-8,'maxiter',100,'tol',dMaxErr,'iErrCompStep',1,'inittype','rand',...
        'itercritupdate',1e-10,'maxorder',1);
    for i=1:iNbUp
        cvsolvFixedRank{i+1} = setparam(cvsolvFixedRank{1},'update',i,'maxorder',100);
    end
    
    if (bManufactRef)
        % Construction d'une solution sur la base des polynomes : accessible pour la PGD
        % sans etre trivial comme les fonctions propres du laplacien
        iNbModes = 50;
        dvOpenCoord = dvCoord(2:end-1); % suppression des noeuds au bord, mieux?
        iDeg = 10;
        imPowers = randi(iDeg+1,iNbModes,iDim)-1;
        dvAlpha = ones(1,iNbModes);
        cmdvModes = cell(iNbModes,iDim);
        for i=1:iNbModes
            for j=1:iDim
                cmdvModes{i,j} = dvOpenCoord.*(1-dvOpenCoord).*dvOpenCoord.^imPowers(i,j);
                dvAlpha(i) = dvAlpha(i)*norm(cmdvModes{i,j});
                cmdvModes{i,j} = cmdvModes{i,j}/norm(cmdvModes{i,j});
            end
        end
        cvsepvRefSol{iSolNum} = SEPMATRIX(cmdvModes,dvAlpha);
        cvsepvSrc{iSolNum} = cvsepmSys{iSolNum}*cvsepvRefSol{iSolNum};
        if (strcmp(strErrCompMode,'reference'))
            for i=1:iNbUp+1
                cvsolvFixedRank{i} = setparam(cvsolvFixedRank{i},'reference',...
                    cvsepvRefSol{iSolNum},'metric',cvsepmSys{iSolNum});
            end
        end
    else
        cvsepvSrc{iSolNum} = SEPMATRIX(cvdvSrc);
    end
    
    for i=1:iNbUp+1
        [cvcvsepvFixedErrSol{i}{iSolNum},cvcvstFixedErrRes{i}{iSolNum}] =...
            solve(cvsepmSys{iSolNum},cvsepvSrc{iSolNum},cvsolvFixedRank{i});
        cvdvFixedErrRank{i}(iSolNum) = cvcvsepvFixedErrSol{i}{iSolNum}.m;
    end
end

figure(round(1/dMaxErr));
plot(2:iMaxDim,cvdvFixedErrRank{1}(1:iNbSols),'r');
hold on
for i=2:iNbUp+1
    plot(2:iMaxDim,cvdvFixedErrRank{i}(1:iNbSols),'Marker',cvstrMarkers{i});
    hold on
end
xlabel('Dimension');
ylabel('Rank');
ylim([0,100]);
legend(cvstrLegends);
strTitle = ['GPGD rank at error ',num2str(iUnits),'e-',num2str(iExp),' vs dimension (',...
    strErrCompMode,' error) - ',strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

if bPrint
    bPrint = false;
    strFileName = lower(strrep(['gpgd-rank-err',num2str(iUnits),'e-',num2str(iExp),'-dim2-',...
        num2str(iMaxDim),'-',strErrCompMode,'-err-',strEllipMode,'-',strTangentialMode],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end


%% Influence du nombre d'updates pour le O-Update (cas degeneres des
%% solutions manufacturees de faible rang)

bDegenerate = true;
if (bDegenerate)
    strKindOfDegen = 'degeneigenf';
    solvCurrent = setparam(solvCurrent,'maxorder',10,'tol',1e-7);
    iNbUp = 1;
    % strKindOfDegen = 'degenpolyf';
    % solvCurrent = setparam(solvCurrent,'maxorder',3,'tol',1e-7);
    % iNbUp = 5;
end

figure(iDim);
% iNbUp = min(iDim+1,length(cvstrMarkers)-1);
% solvNUp = setparam(solvCurrent,'update',27);
% [cvsepvNUpSols{6},cvstNUpResults{6}] = solve(cvsepmSys,cvsepvSrc,solvNUp);
% bvNUpSolsComputed(6) = true;
% semilogy(cvstNUpResults{6}.error,'Marker',cvstrMarkers{6});
% cvstrLegends{6} = ['(',num2str(27),')G-PGD'];
% bvNUpSolsComputed(1:iNbUp+1)=0;
if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,iNbUp+1);
    cvsepvNUpSols = cell(1,iNbUp);
    cvstNUpResults = cell(1,iNbUp);
end
for i=0:iNbUp
    if (~bvNUpSolsComputed(i+1))
        solvNUp = setparam(solvCurrent,'update',i);
        solvNUp = setparam(solvCurrent,'update',27);
        [cvsepvNUpSols{i+1},cvstNUpResults{i+1}] = solve(cvsepmSys,cvsepvSrc,solvNUp);
        bvNUpSolsComputed(i+1) = true;
    end
    if (~i)
        semilogy(cvstNUpResults{i+1}.error,'r','Marker',cvstrMarkers{i+1});
    else
        semilogy(cvstNUpResults{i+1}.error,'Marker',cvstrMarkers{i+1});
    end
    hold on;
    if (~i)
        % cvstrLegends{i+1} = 'G-PGD';
        if (~bvNSTUpSolsComputed(1))
            cvsepvNSTUpSols{1} = cvsepvNUpSols{i+1};
            cvstNSTUpResults{1} = cvstNUpResults{i+1};
            bvNSTUpSolsComputed(1) = true;
        end
    elseif (i == iNbUp)
        % cvstrLegends{i+1} = ['(',num2str(i),')G-PGD'];
        if (~bvNSTUpSolsComputed(2))
            cvsepvNSTUpSols{2} = cvsepvNUpSols{i+1};
            cvstNSTUpResults{2} = cvstNUpResults{i+1};
            bvNSTUpSolsComputed(2) = true;
        end
    else
        % cvstrLegends{i+1} = ['(',num2str(i),')G-PGD'];
    end
end

xlabel('Rank');
ylabel('Error');
legend(cvstrLegends);

if (bDegenerate)
    strTitle = [strDim,'GPGD convergence (',strErrCompMode,' error) - ',...
        strEllipMode,' ',strTangentialMode,strKindOfDegen,' problem'];
    title(strTitle);
    strFileName = lower(strrep([strDim,'nup-GPGD-cv-',strErrCompMode,'-err-',...
        strEllipMode,'-',strTangentialMode,'-',strKindOfDegen],' ','-'));
else
    strTitle = [strDim,'GPGD convergence (',strErrCompMode,' error) - ',...
        strEllipMode,' ',strTangentialMode,' problem'];
    title(strTitle);
    strFileName = lower(strrep([strDim,'nup-GPGD-cv-',strErrCompMode,'-err-',...
        strEllipMode,'-',strTangentialMode],' ','-'));
end

if bPrint
    bPrint = false;
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end


%% Evolution des modes vs updates
bRes = false;
bAdj = false;
iMaxOrder = iNbModes;
for i=1:iNbModes*iDim
    figure(i);
end
for i=0:20
    solvTest = setparam(solvCurrent,'maxorder',iMaxOrder,'itercrit',5e-8,...
        'maxiter',100,'tol',1e-7,'update',0,'inittype','rand',...
        'errorindicator','reference','reference',cvsepvRefSol,'itercritupdate',1e-10,...
        'iNonLinearUpdates',0,'ivModesToUpdate',1:iMaxOrder,'residual',bRes,'adjoint',bAdj,...
        'update',i);
    [sepvTestSol,stResult] = solve(cvsepmSys,cvsepvSrc,solvTest);
    for j=1:iNbModes
        for k=1:iDim
            set(0,'CurrentFigure',(j-1)*iNbModes+k);
            % plot(dvOpenCoord,sepvTestSol.F{j,k});
            hold on;
        end
    end
end
% figure(1)
% if (bRes)
%     semilogy(stResult.error,'r');
% elseif (bAdj)
%     semilogy(stResult.error,'g');
% else
%     semilogy(stResult.error);
% end
% hold on
% xlabel('Rank');
% ylabel('Error');
% strTitle = [strDim,'GPGD convergence (',strErrCompMode,' error) - ',...
%     strEllipMode,' ',strTangentialMode,' problem'];
% title(strTitle);


%% Convergence au rang d'entree
bRes = false;
bAdj = false;
iMaxOrder = iNbModes;
iMaxUp = 1000;
solvTest = setparam(solvCurrent,'maxorder',iMaxOrder,'itercrit',5e-8,...
    'maxiter',100,'tol',1e-7,'update',0,'inittype','rand',...
    'errorindicator','reference','reference',cvsepvRefSol,'itercritupdate',1e-10,...
    'iNonLinearUpdates',0,'ivModesToUpdate',1:iMaxOrder,'residual',bRes,'adjoint',bAdj,...
    'update',iMaxUp,'updatestep',iMaxOrder);
[sepvTestSol,stResult] = solve(cvsepmSys,cvsepvSrc,solvTest);
figure(iMaxUp)
if (bRes)
    semilogy(stResult.error,'r');
elseif (bAdj)
    semilogy(stResult.error,'g');
else
    semilogy(stResult.error);
end
hold on
xlabel('Rank');
ylabel('Error');
strTitle = [strDim,'GPGD convergence (',strErrCompMode,' error) - ',...
    strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

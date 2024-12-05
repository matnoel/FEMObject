% Applications numeriques pour la publication "Tensor-based methods and
% Proper Generalized Decompositions for the numerical solution of high
% dimensional problems: alternative definitions", G. Bonithon, A. Nouy,
% 2011.

% Exemple 2 :
% Equation d'advection-diffusion -\Delta u + <\beta(x_1,...,x_n),\nabla u> = f
% Analyse du comportement de la G-PGD et de l'impact des updates par rapport
% a l'elevation en dimension, l'ellipticite, et la tangentialite de
% l'operateur

%% Preambule

clear();
bPrint = false;
bElliptic = false;
bTangential = true;
bOne = false;
bInvertNonSym = false;
bvNUpSolsComputed = zeros(1,1000);
strErrCompMode = 'residual';
dFreq = 10;
dAmp = 30;
iExp = 0;
dDiff = 10^(-iExp);
iDim = 3;
iNbElem = 1000;
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
strPath = [pwd,'/Publi/Ex2/'];

% Generation du modele
domD1d = DOMAIN(1);
modM1d = mesh(domD1d,iNbElem);
modM1d = createddlnode(modM1d,DDL('u'));
modM1d = addcl(modM1d,POINT([0;1]),'u',0);

% Terme de diffusion (Id\otimes\Delta) et terme source pour une dimension
% Formulation continue
dvCoord = getcoord(getnode(modM1d));
blinfD00 = BILINFORM(0,0);
dvDiff = dDiff*ones(size(dvCoord,1),1);
blinfD11 = BILINFORM(1,1,dvDiff,0);
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

% Terme d'advection
% Formulation continue
if (bOne)
    cmdvWeight(1:iDim,1:iDim) = {-ones(size(dvCoord))};
else
    cmdvWeight(1:iDim,1:iDim) = {dAmp*cos(dFreq*pi*dvCoord)};
    % cmdvWeight(1:iDim,1:iDim) = {cos(pi*dvCoord)};
    % cmdvWeight(1:iDim,1:iDim) = {dvCoord ones(size(dvCoord));...
    %     ones(size(dvCoord)) dvCoord};
end

if (bElliptic)
    cmdvWeight = cellfun(@(x)abs(x),cmdvWeight,'UniformOutput',false);
end
cmblinfD00 = cellfun(@(x)BILINFORM(0,0,x,0),cmdvWeight,'UniformOutput',false);
cmblinfD10 = cellfun(@(x)BILINFORM(1,0,x,0),cmdvWeight,'UniformOutput',false);
cmblinfD01 = cellfun(@(x)BILINFORM(0,1,x,0),cmdvWeight,'UniformOutput',false);
if (bTangential)
    blinfD00one = BILINFORM(0,0);
end

% Formulation discrete
stvS(1).type = '{}';
stvS(1).subs = {modM1d};
stvS(2).type = '()';
stvS(2).subs = {':',':'};
cmdmD00 = cell(iDim,iDim);
cmdmD10 = cell(iDim,iDim);
cmdmD01 = cell(iDim,iDim);
for i=1:iDim
    for j=1:iDim
        cmdmD00{i,j} = cmblinfD00{i,j}{modM1d}(:,:);
        cmdmD10{i,j} = cmblinfD10{i,j}{modM1d}(:,:);
        cmdmD01{i,j} = cmblinfD01{i,j}{modM1d}(:,:);
    end
end
if (bTangential)
    dmD00one = blinfD00one{modM1d}(:,:);
end

% Assemblage des dimensions
if (bTangential)
    cmdmSys(iDim+1:2*iDim,1:iDim) = {dmD00one};
else
    cmdmSys(iDim+1:2*iDim,1:iDim) = cmdmD00;
end
if (bInvertNonSym)
    cmdmSys(logical([zeros(iDim);eye(iDim)])) = cmdmD10(logical(eye(iDim)));
else
    cmdmSys(logical([zeros(iDim);eye(iDim)])) = cmdmD01(logical(eye(iDim)));
end

% Passage aux sepmatrix
sepmSys = SEPMATRIX(cmdmSys);
sepvSrc = SEPMATRIX(cvdvSrc);

% Solveur par defaut
solvDefault = SEPSOLVER(iDim,'adjoint',0,'inittype','one',...
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

solvCurrent = setparam(solvCurrent,'maxorder',10,'inittype','rand');
iNbUp = 1;
bAdj = 0;
if (bAdj)
    strPrePGD = 'MM';
else
    strPrePGD = 'G';
end
figure(iDim+bAdj);

if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,iNbUp+1);
    cvsepvNUpSols = cell(1,iNbUp);
    cvstNUpResults = cell(1,iNbUp);
else
    % bvNUpSolsComputed(1:end)=0;
end
for i=0:iNbUp
    if (~bvNUpSolsComputed(i+1))
        solvNUp = setparam(solvCurrent,'update',i,'adjoint',bAdj);
        [cvsepvNUpSols{i+1},cvstNUpResults{i+1}] = solve(sepmSys,sepvSrc,solvNUp);
        bvNUpSolsComputed(i+1) = true;
    end
    if (~i)
        semilogy(cvstNUpResults{i+1}.error,'r','Marker',cvstrMarkers{i+1});
        cvstrLegends{i+1} = [strPrePGD,'-PGD'];
    else
        semilogy(cvstNUpResults{i+1}.error,'Marker',cvstrMarkers{i+1});
        cvstrLegends{i+1} = ['(',num2str(i),')',strPrePGD,'-PGD'];
    end
    hold on;
end

xlabel('Rank');
ylabel('Error');
legend(cvstrLegends,'Location','SouthWest');

strTitle = [strDim,strPrePGD,'-PGD convergence (',strErrCompMode,' error) - ',...
    strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,'nup-',strPrePGD,'-PGD-cv-',strErrCompMode,'-err-',...
        strEllipMode,'-',strTangentialMode],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end

% figure;
% temp = expand(cvsepvNUpSols{end})';
% [S,X,Y ] = mesh(DOMAIN(2),iNbElem,iNbElem);
% S = createddlnode(S,DDL('u'));
% S = addcl(S,[],'u');
% plot(temp(:),S)%,'surface')
% colorbar;
% title('Solution, separated form')


%% Diffusivity parameter influence

bTangential = true;
bElliptic = true;
solvCurrent = setparam(solvCurrent,'maxorder',20);
iNbUp = 1;
bAdj = 1;
if (bAdj)
    strPrePGD = 'MM';
else
    strPrePGD = 'G';
end

figure(iDim+bAdj);
iMaxExp = 4;
if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,iMaxExp+1);
    cvsepvNUpSols = cell(1,(iMaxExp+1)*(iNbUp+1));
    cvstNUpResults = cell(1,(iMaxExp+1)*(iNbUp+1));
else
    % bvNUpSolsComputed(1:end)=0;
end
iNbSols = 0;

for l=0:iNbUp
    for k=0:iMaxExp
        iNbSols = iNbSols + 1;
        if (~k)
            cvstrLegends{k+1} = '\epsilon = 1';
        else
            cvstrLegends{k+1} = ['\epsilon = 1e-',num2str(k)];
        end
        % Terme de diffusion (Id\otimes\Delta) et terme source pour une dimension
        % Formulation continue
        dDiff = 10^(-k);
        dvCoord = getcoord(getnode(modM1d));
        blinfD00 = BILINFORM(0,0);
        dvDiff = dDiff*ones(size(dvCoord,1),1);
        blinfD11 = BILINFORM(1,1,dvDiff,0);
        dvSrc = ones(size(dvCoord,1),1);
        linfSrc = LINFORM(0,dvSrc,0);
        
        % Formulation discrete
        dmD00 = blinfD00{modM1d}(:,:);
        dmD11 = blinfD11{modM1d}(:,:);
        dvSrc1d = linfSrc{modM1d}{:};
        
        % Assemblage des dimensions
        cvdvSrc(1:iDim) = {dvSrc1d};
        cmdmSys = cell(iDim,iDim);
        cmdmSys(1:iDim,1:iDim) = {dmD00};
        cmdmSys(logical(eye(iDim))) = {dmD11};
        
        % Terme d'advection
        % Formulation continue
        cmdvWeight(1:iDim,1:iDim) = {ones(size(dvCoord))};
        if (bElliptic)
            cmdvWeight = cellfun(@(x)abs(x),cmdvWeight,'UniformOutput',false);
        end
        cmblinfD00 = cellfun(@(x)BILINFORM(0,0,x,0),cmdvWeight,'UniformOutput',false);
        cmblinfD10 = cellfun(@(x)BILINFORM(1,0,x,0),cmdvWeight,'UniformOutput',false);
        if (bTangential)
            blinfD00one = BILINFORM(0,0);
        end
        
        % Formulation discrete
        stvS(1).type = '{}';
        stvS(1).subs = {modM1d};
        stvS(2).type = '()';
        stvS(2).subs = {':',':'};
        cmdmD00 = cell(iDim,iDim);
        cmdmD10 = cell(iDim,iDim);
        for i=1:iDim
            for j=1:iDim
                cmdmD00{i,j} = cmblinfD00{i,j}{modM1d}(:,:);
                cmdmD10{i,j} = cmblinfD10{i,j}{modM1d}(:,:);
            end
        end
        if (bTangential)
            dmD00one = blinfD00one{modM1d}(:,:);
        end
        
        % Assemblage des dimensions
        if (bTangential)
            cmdmSys(iDim+1:2*iDim,1:iDim) = {dmD00one};
        else
            cmdmSys(iDim+1:2*iDim,1:iDim) = cmdmD00;
        end
        cmdmSys(logical([zeros(iDim);eye(iDim)])) = cmdmD10(logical(eye(iDim)));
        
        % Passage aux sepmatrix
        sepmSys = SEPMATRIX(cmdmSys);
        sepvSrc = SEPMATRIX(cvdvSrc);
        
        if (~bvNUpSolsComputed(iNbSols))
            solvNUp = setparam(solvCurrent,'update',l,'adjoint',bAdj);
            [cvsepvNUpSols{iNbSols},cvstNUpResults{iNbSols}] = solve(sepmSys,sepvSrc,solvNUp);
            bvNUpSolsComputed(iNbSols) = true;
        end
        if (~l)
            semilogy(cvstNUpResults{iNbSols}.error,'r','Marker',cvstrMarkers{k+1});
        else
            semilogy(cvstNUpResults{iNbSols}.error,'--','Marker',cvstrMarkers{k+1});
        end
        hold on;
    end
end

xlabel('Rank');
ylabel('Error');
legend(cvstrLegends,'Location','SouthWest');

strTitle = [strDim,strPrePGD,'-PGD convergence (',strErrCompMode,' error) - ',...
    strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,'nup-',strPrePGD,'-PGD-cv-',strErrCompMode,'-err-',...
        strEllipMode,'-',strTangentialMode],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end


%% Superposition de G-PGD et MM-PGD a \epsilon fixe

solvCurrent = setparam(solvCurrent,'maxorder',20);
iNbUp = 1;
bAdj = 1;
figure(iDim+bAdj);

if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,(iNbUp+1)*(bAdj+1));
    cvsepvNUpSols = cell(1,(iNbUp+1)*(bAdj+1));
    cvstNUpResults = cell(1,(iNbUp+1)*(bAdj+1));
else
    %bvNUpSolsComputed(1:end)=0;
end
iNbSols = 0;
for i=0:iNbUp
    for j=0:bAdj
        iNbSols = iNbSols + 1;
        if (~bvNUpSolsComputed(iNbSols))
            solvNUp = setparam(solvCurrent,'update',i,'adjoint',j);
            [cvsepvNUpSols{iNbSols},cvstNUpResults{iNbSols}] = solve(sepmSys,sepvSrc,solvNUp);
            bvNUpSolsComputed(iNbSols) = true;
        end
        if (~i)
            if (~j)
                semilogy(cvstNUpResults{iNbSols}.error,'r','Marker',cvstrMarkers{i+1});
                strPrePGD = 'G';
            else
                semilogy(cvstNUpResults{iNbSols}.error,'--r','Marker',cvstrMarkers{i+1});
                strPrePGD = 'MM';
            end
            cvstrLegends{iNbSols} = [strPrePGD,'-PGD'];
        else
            if (~j)
                semilogy(cvstNUpResults{iNbSols}.error,'Marker',cvstrMarkers{i+1});
                strPrePGD = 'G';
            else
                semilogy(cvstNUpResults{iNbSols}.error,'Marker',cvstrMarkers{i+1},'Linestyle','--');
                strPrePGD = 'MM';
            end
            cvstrLegends{iNbSols} = ['(',num2str(i),')',strPrePGD,'-PGD'];
        end
        hold on;
    end
end

xlabel('Rank');
ylabel('Error');
legend(cvstrLegends,'Location','SouthWest');

strTitle = [strDim,'G-MM-PGD convergence at \epsilon = 1e-',num2str(iExp),...
    ' (',strErrCompMode,' error) - ',strEllipMode,' ',strTangentialMode,' problem'];
title(strTitle);

if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,'nup-','G-MM-PGD-cv-eps-1e-',num2str(iExp),...
        strErrCompMode,'-err-',strEllipMode,'-',strTangentialMode],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end

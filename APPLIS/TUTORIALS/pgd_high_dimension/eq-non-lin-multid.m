% Applications numeriques pour la publication "Tensor-based methods and
% Proper Generalized Decompositions for the numerical solution of high 
% dimensional problems: alternative definitions", G. Bonithon, A. Nouy,
% 2011.

% Exemple 3 et 4 :
% Equation de diffusion reaction non lineaire -\epsilon\Delta u + \alpha(x_1,...,x_d)u^3 = f
% Comparaison des approches Newton-G-PGD et G-PGD-Newton
% Equation d'advection-diffusion-reaction non lineaire
% -\epsilon\Delta u + div(u^n\nabla u) = f
% Test uniquement de l'approche Newton-PGD mais en G-PGD et MM-PGD

%% Preambule
clear();
bPrint = false;
bIncrement = true;
bvNUpSolsComputed = zeros(1,1000);
strErrCompMode = 'residual';
dFreq = 10;
dAmp = 30;
iExp = 1;
dDiff = 10^(-iExp);
iDim = 2;
iNbElem = 1000;
strDim = [num2str(iDim),'D-'];
strPath = [pwd,'/Publi/Ex3/'];

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

iNewtPow = 2;
bInitialize = false;
cmdvSrc = {-sin(pi*dvCoord).^(iNewtPow+1) pi^2*sin(pi*dvCoord).^(iNewtPow-1).*...
    (iNewtPow.*cos(pi*dvCoord).^2-sin(pi*dvCoord).^2);...
    pi^2*sin(pi*dvCoord).^(iNewtPow-1).*...
    (iNewtPow.*cos(pi*dvCoord).^2-sin(pi*dvCoord).^2) ...
    -sin(pi*dvCoord).^(iNewtPow+1)};
cmlinfSrc = cellfun(@(x)LINFORM(0,x,0),cmdvSrc,'UniformOutput',false);

% Formulation discrete
dmD00 = blinfD00{modM1d}(:,:);
dmD11 = blinfD11{modM1d}(:,:);
dvSrc1d = linfSrc{modM1d}{:};
for i=1:2
    for j=1:2
        cmdvSrc{i,j} = cmlinfSrc{i,j}{modM1d}{:};
    end
end

% Assemblage des dimensions
cvdvSrc(1:iDim) = {dvSrc1d};
cmdmSys(1:iDim,1:iDim) = {dmD00};
cmdmSys(logical(eye(iDim))) = {dmD11};

% Terme de reaction, formulation continue
blinfD00 = BILINFORM(0,0);

% Formulation discrete
dmD00 = blinfD00{modM1d}(:,:);

% Assemblage des dimensions
cmdmSys(iDim+1,1:iDim) = {dmD00};

% Passage aux sepmatrix
sepmSys = SEPMATRIX(cmdmSys);

if (bInitialize)
    sepvSrc = SEPMATRIX(cmdvSrc);
else
    sepvSrc = SEPMATRIX(cvdvSrc);
end

% Solveur par defaut
solvDefault = SEPSOLVER(iDim,'adjoint',0,'inittype','one',...
    'residual',0,'maxorder',20,'tol',1e-6,'update',0,....
    'errorindicator',strErrCompMode,'ortho',0,'updateeps',0,...
    'itercrit',5e-3);
cvstrMarkers = {'none','+','*','o','^','s','p','x',...
    'none','+','*','o','^','s','p','x'};

solvCurrent = setparam(solvDefault,'maxorder',12,'itercrit',5e-8,...
    'maxiter',100,'tol',1e-7,'iErrCompStep',1,'inittype','rand',...
    'itercritupdate',1e-10);


%% Newton non symetrique (exemple 4)

iNbUp = 0;
bAdj = 0;
iNbIterNewton = 1;
iMaxRank = 10;
iNbUpSvd = 2;
iRankSvd = 3;
dErrSvd = 1e-6;
if (bAdj)
    strPrePGD = 'MM';
else
    strPrePGD = 'G';
end
if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,iNbIterNewton);
    cvsepvNUpSols = cell(iMaxRank,iNbIterNewton);
    cvstNUpResults = cell(iMaxRank,iNbIterNewton);
    dvNewtonErr = zeros(iMaxRank,iNbIterNewton);
    sepvModifSrc = sepvSrc;
    sepmInitSys = sepmSys;
    sepvInitSrc = sepvSrc;
else
%bvNUpSolsComputed(1:end,1:end)=0;
end

dRefNorm = norm(sepvSrc);
figure(iNbUp+bAdj+1);
for r=1:iMaxRank
    if (~bvNUpSolsComputed(r))
        bvNUpSolsComputed(r) = true;
        solvCurrent = setparam(solvCurrent,'maxorder',r,'inittype','rand','update',iNbUp,'adjoint',bAdj,...
            'bNlSolver',0);
        for i=1:iNbIterNewton
% Initialisation forcee
            if (bInitialize && i==1)
                cmdvSol{1} = sin(pi*(dvCoord(2:end-1)+0.1));
%cmdvSol{1} = ones(size(dvCoord(2:end-1)));
%cmdvSol{1} = cmdvSol{1}/norm(cmdvSol{1});
                cmdvSol{2} = ones(size(dvCoord(2:end-1)));
                cmdvSol(2) = cmdvSol(1);
%     cmdvSol(1) = cmdvSol(2);
                cvsepvNUpSols{r,i} = SEPMATRIX(cmdvSol);
%cvsepvNUpSols{i} = seprand((iNbElem-1)*ones(size(dvCoord)));
                cvstNUpResults{r,i}.error = zeros(1,r);
            elseif (i==1)
                disp('a')
                [cvsepvNUpSols{r,i},cvstNUpResults{r,i}] = solve(sepmInitSys,sepvInitSrc,solvCurrent);
            else
% Calcul de la solution approchee
                [cvsepvNUpSols{r,i},cvstNUpResults{r,i}] = solve(sepmSys,sepvModifSrc,solvCurrent);
            end

            if (bIncrement)
                if (i==1)
                    sepvPastSol = cvsepvNUpSols{r,i};
                else
                    sepvPastSol = cvsepvNUpSols{r,i} + sepvPastSol;
                end
            else
                sepvPastSol = cvsepvNUpSols{r,i};
            end

            if (iRankSvd<r*i)
                solvSvd = setparam(solvCurrent,'maxorder',iRankSvd,'update',iNbUpSvd,...
                    'errorindicator','residual','tol',dErrSvd,'display',false);
                [sepvRedPastSol,stRes] = multisvd(sepvPastSol,solvSvd);
                dvSvdErr(r,i) = stRes.error(end);
            else
                sepvRedPastSol = sepvPastSol;
                dvSvdErr(r,i) = 0;
            end

% Actualisation du systeme
% Formulation continue
            cmdvWeight = sepvRedPastSol.F;
            sepvSolPow = sepvRedPastSol;
            if (iNewtPow == 1)
                cmdvSolPowM1(1,1:iDim) = {ones(size(dvCoord))};
            else
                cmdvSolPowM1(1,1:iDim) = {ones(size(dvCoord(2:end-1)))};
            end
            sepvSolPowM1 = SEPMATRIX(cmdvSolPowM1);
            for j=2:iNewtPow
                sepvSolPow = sepvSolPow.*sepvRedPastSol;
                sepvSolPowM1 = sepvSolPowM1.*sepvRedPastSol;
            end
            iNbModes = sepvSolPow.m;
            iNbModes0 = sepvRedPastSol.m;
            iNbModes1 = sepvSolPowM1.m;
            cmdvWeightPow = sepvSolPow.F;
            cmdvWeightPowM1 = sepvSolPowM1.F;
            cmblinfD00Lap = BILINFORM(0,0);
            cmblinfD11Lap = BILINFORM(1,1,dvDiff,0);
            cmblinfD00 = cellfun(@(x)BILINFORM(0,0,x,0),cmdvWeightPow,'UniformOutput',false);
            cmblinfD11 = cellfun(@(x)BILINFORM(1,1,x,0),cmdvWeightPow,'UniformOutput',false);
            cmblinfD01 = cellfun(@(x)MULTILINFORM([1,1,0,0]),cmdvWeightPow,'UniformOutput',false);

% Formulation discrete
            cmdmSymSys = cell(iNbModes*iDim,iDim);
            cmdmNSymSys = cell(iNbModes*iDim,iDim);
            cmdmLapSys = cell(iDim,iDim);
            dvAlpha = ones(1,iNbModes*iDim);
            for mm=1:iNbModes0
                for mmm=1:iNbModes1
                    m = (mm-1)*iNbModes1+mmm;
                    dvAlpha(1,1+(m-1)*iDim:m*iDim) = sepvSolPowM1.alpha(mmm)*sepvRedPastSol.alpha(mm);
                    for d=1:iDim
                        for dd=1:iDim
                            if (d==dd)
                                cmdmLapSys{dd,d} = cmblinfD11Lap{modM1d}(:,:);
                                cmdmSymSys{dd+(m-1)*iDim,d} = cmblinfD11{m,d}{modM1d}(:,:);
                                cmdmNSymSys{dd+(m-1)*iDim,d} = cmblinfD01{m,d}{modM1d}(:,cmdvWeight{mm,d},...
                                    cmdvWeightPowM1{mmm,d},:);
                            else
                                cmdmLapSys{dd,d} = cmblinfD00Lap{modM1d}(:,:);
                                cmdmSymSys{dd+(m-1)*iDim,d} = cmblinfD00{m,d}{modM1d}(:,:);
                                cmdmNSymSys{dd+(m-1)*iDim,d} = cmblinfD00{m,d}{modM1d}(:,:);
                            end
                        end
                    end
                end
            end

% Passage aux sepmatrix
            sepmLapSys = SEPMATRIX(cmdmLapSys);
            sepmSymSys = SEPMATRIX(cmdmSymSys,dvAlpha);
            sepmNSymSys = SEPMATRIX(cmdmNSymSys,iNewtPow*dvAlpha);
            sepmSys = sepmLapSys + sepmSymSys + sepmNSymSys;
            sepvPerturbSrc = sepmNSymSys*sepvRedPastSol;
            sepvNonLinSol = (sepmLapSys+sepmSymSys)*sepvPastSol;
            if (bIncrement)
                sepvModifSrc = sepvSrc - sepvNonLinSol;
            else
                sepvModifSrc = sepvSrc + sepvPerturbSrc;
            end

% Calcul de l'erreur pour l'algo de Newton
            dvNewtonErr(r,i) = norm(sepvSrc-sepvNonLinSol)/dRefNorm;
        end
    end

% Trace de l'erreur pour l'algo de Newton
    cvstrLegends{r} = ['m = ',num2str(r)];
    if (~iNbUp)
        semilogy(dvNewtonErr(r,1:iNbIterNewton),'r','Marker',cvstrMarkers{r});
    else
        semilogy(dvNewtonErr(r,1:iNbIterNewton),'Marker',cvstrMarkers{r});
    end
    hold on
end

xlabel('k');
ylabel('Error');
legend(cvstrLegends(1:iMaxRank),'Location','SouthWest');
strTitle = [strDim,'Newton-convergence (',strErrCompMode,' error) - ',...
    'non elliptic problem'];
title(strTitle);
if bPrint
    strFileName = lower(strrep([strDim,'Newton-cv-',strErrCompMode,'-err-',...
        num2str(iNbUp),strPrePGD,'-PGD-non-elliptic'],' ','-'));
    title('');
    myprint(strPath,strFileName,'epsc2');
end

% Trace des courbes PGD
figure(iNbUp+bAdj+10);
for i=1:iNbIterNewton
    semilogy(cvstNUpResults{iMaxRank,i}.error,'Marker',cvstrMarkers{i});
    if (iNbUp)
        cvstrLegends{i} = ['(',num2str(iNbUp),')',strPrePGD,'-PGD, k=',num2str(i)];
    else
        cvstrLegends{i} = [strPrePGD,'-PGD, k=',num2str(i)];
    end
    hold on;
end
xlabel('Rank');
ylabel('Error');
legend(cvstrLegends(1:iNbIterNewton),'Location','SouthWest');
strTitle = [strDim,'Newton-',num2str(iNbUp),strPrePGD,'-PGD convergence (',strErrCompMode,' error) - ',...
    'non elliptic problem'];
title(strTitle);
if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,'Newton-',num2str(iNbUp),strPrePGD,'-PGD-cv-',...
        strErrCompMode,'-err-','non-elliptic'],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end


%% PGD symetrique (exemple 3)

iNbUp = 0;
bAdj = 0;
iMaxRank = 5;
if (bAdj)
    strPrePGD = 'MM';
else
    strPrePGD = 'G';
end
%bForce = true;
if (~exist('sepvSol','var') || bForce)
    solvCurrent = setparam(solvCurrent,'maxorder',iMaxRank,'inittype','one','update',iNbUp,'adjoint',bAdj,...
        'bNlSolver',1e-6,'cmdmLinOp',cmdmSys(1:iDim,1:iDim),'modM1d',modM1d,'maxiter',100,'itercrit',5e-3);
    [sepvSol,stResult] = solve(sepmSys,sepvSrc,solvCurrent);
    bForce = false;
end
semilogy(stResult.error,'b');
hold on;
strFileName = [strPath,num2str(iDim),'d-newton-cv-residual-err-',num2str(iNbUp),'g-pgd-elliptic.mat'];
load(strFileName,'cvstNUpResults');
semilogy(cvstNUpResults{end}.error);
xlim([1,iMaxRank]);
if (iNbUp)
    cvstrLegends{1} = ['Progressive (',num2str(iNbUp),')',strPrePGD,'-PGD'];
    cvstrLegends{2} = ['Newton+(',num2str(iNbUp),')',strPrePGD,'-PGD'];
else
    cvstrLegends{1} = ['Progressive ',strPrePGD,'-PGD'];
    cvstrLegends{2} = ['Newton+',strPrePGD,'-PGD'];
end
legend(cvstrLegends);

xlabel('Rank');
ylabel('Error');
strTitle = [strDim,num2str(iNbUp),strPrePGD,'-PGD-Newton convergence (',strErrCompMode,' error) - ',...
    'elliptic problem'];
title(strTitle);
if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,num2str(iNbUp),strPrePGD,'-PGD-Newton-cv-',strErrCompMode,'-err-',...
        'elliptic'],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end



%% Newton symetrique (exemple 3)

iNbUp = 0;
bAdj = 0;
iNbIterNewton = 5;
iMaxRank = 10;
if (bAdj)
    strPrePGD = 'MM';
else
    strPrePGD = 'G';
end
if (~exist('cvsepvNUpSols','var'))
    cvstrLegends = cell(1,iNbIterNewton);
    cvsepvNUpSols = cell(1,iNbIterNewton);
    cvstNUpResults = cell(1,iNbIterNewton);
    dvNewtonErr = zeros(1,iNbIterNewton);
    sepvModifSrc = sepvSrc;
else
%bvNUpSolsComputed(1:end)=0;
end

iNbUpSvd = 2;
iRankSvd = iMaxRank;
dErrSvd = 1e-6;
dRefNorm = norm(sepvSrc);
solvCurrent = setparam(solvCurrent,'maxorder',iMaxRank,'inittype','rand','update',iNbUp,'adjoint',bAdj,...
    'bNlSolver',0);
figure(iNbUp+bAdj+1);
for i=1:iNbIterNewton
    if (~bvNUpSolsComputed(i))
% Calcul de la solution approchee
        [cvsepvNUpSols{i},cvstNUpResults{i}] = solve(sepmSys,sepvModifSrc,solvCurrent);
        bvNUpSolsComputed(i) = true;

        if (bIncrement)
            if (i==1)
                sepvPastSol = cvsepvNUpSols{i};
            else
                sepvPastSol = cvsepvNUpSols{i} + sepvPastSol;
            end
        else
            sepvPastSol = cvsepvNUpSols{i};
        end

% Actualisation du systeme, terme de reaction, formulation continue
        sepvUWeight = sepvPastSol.*sepvPastSol;
        solvSvd = setparam(solvCurrent,'maxorder',iRankSvd,'update',iNbUpSvd,...
            'errorindicator','residual','tol',dErrSvd);
        sepvUWeight = multisvd(sepvUWeight,solvSvd);
        cmdvWeight = sepvUWeight.F;
        cmblinfD00 = cellfun(@(x)BILINFORM(0,0,x,0),cmdvWeight,'UniformOutput',false);

% Formulation discrete
        cmdmD00 = cell(size(cmdvWeight));
        for m=1:size(cmdvWeight,1)
            for d=1:iDim
                cmdmD00{m,d} = cmblinfD00{m,d}{modM1d}(:,:);
            end
        end

% Assemblage des dimensions
        cmdmSys = [cmdmSys(1:iDim,1:iDim);cell(size(cmdvWeight,1),iDim)];
        cmdmSys(iDim+1:iDim+size(cmdvWeight,1),1:iDim) = cmdmD00;

% Passage aux sepmatrix
        dvAlpha = [ones(1,iDim),3*sepvUWeight.alpha]; % derivation de l'operateur
        sepmSys = SEPMATRIX(cmdmSys,dvAlpha);
        dvAlpha = 2*sepvUWeight.alpha; % derivation de l'operateur - 1
        sepmPerturbSys = SEPMATRIX(cmdmD00,dvAlpha);
        sepvNonLinSol = sepmSys*sepvPastSol - sepmPerturbSys*sepvPastSol;
        if (bIncrement)
            sepvModifSrc = sepvSrc - sepvNonLinSol;
        else
            sepvModifSrc = sepvSrc + sepmPerturbSys*sepvPastSol;
        end
%         solvSvd = setparam(solvCurrent,'maxorder',iRankSvd,'update',iNbUpSvd,...
%             'errorindicator','residual','tol',dErrSvd);
%         sepvModifSrc = multisvd(sepvModifSrc,solvSvd);

% Calcul de l'erreur pour l'algo de Newton
        dvNewtonErr(i) = norm(sepvSrc-sepvNonLinSol)/dRefNorm;
    end

% Trace des courbes
    semilogy(cvstNUpResults{i}.error,'Marker',cvstrMarkers{i});
    if (iNbUp)
        cvstrLegends{i} = ['(',num2str(iNbUp),')',strPrePGD,'-PGD, k=',num2str(i)];
    else
        cvstrLegends{i} = [strPrePGD,'-PGD, k=',num2str(i)];
    end
    hold on;
end

xlabel('Rank');
ylabel('Error');
legend(cvstrLegends);
strTitle = [strDim,'Newton-',num2str(iNbUp),strPrePGD,'-PGD convergence (',strErrCompMode,' error) - ',...
    'elliptic problem'];
title(strTitle);
if bPrint
    strFileName = lower(strrep([strDim,'Newton-',num2str(iNbUp),strPrePGD,'-PGD-cv-',strErrCompMode,'-err-',...
        'elliptic'],' ','-'));
    title('');
    myprint(strPath,strFileName,'epsc2');
end

% Trace de l'erreur pour l'algo de Newton
figure(iNbUp+bAdj+10);
semilogy(dvNewtonErr(1:iNbIterNewton),'r')
xlabel('k');
ylabel('Error');
if (iDim==5 && iNbUp==1)
    ylim([5e-4,1]);
end
strTitle = [strDim,'Newton-convergence (',strErrCompMode,' error) - ',...
    'elliptic problem'];
title(strTitle);
if bPrint
    bPrint = false;
    strFileName = lower(strrep([strDim,'Newton-cv-',strErrCompMode,'-err-',...
        num2str(iNbUp),strPrePGD,'-PGD-elliptic'],' ','-'));
    title('');
    save([strPath,strFileName]);
    myprint(strPath,strFileName,'epsc2');
end

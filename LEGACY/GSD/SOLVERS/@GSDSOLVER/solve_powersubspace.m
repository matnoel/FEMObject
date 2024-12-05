function [u,result] = solve_powersubspace(GSD,A,b,ureuse,varargin)
% function u = solve_powersubspace(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

display_ = getparam(GSD,'display');

fprintf('\n GSD solveur ... POWER-SUBSPACE')
if display_
    fprintf('\n')
end
paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
fprintf('DIMENSION M = %d\n',paramradial.nbfoncmax)

paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
paramradial.update = getparam(GSD,'update');
orthocrit = getparam(GSD,'orthocrit');

errorindicator = getparam(GSD,'errorindicator');
uref = getcharin('reference',varargin);

if ~isempty(uref)
    errorindicator='reference';
end

[A,b,PC]=pcsystemupdate(A,b,varargin{:});

if getparam(GSD,'righthandSD');
    tolfactSD = getparam(GSD,'finalSDfacttol');
    
    if isa(b,'PCRADIALMATRIX')
        mm=length(b);
    else
        mm=100;
    end
    b = spectral_decomposition(b,'tol',paramradial.tol*tolfactSD,...
        'nbfoncmax',mm);
end

P = getP(PC);
n=size(A,1);

errorpf=zeros(1,paramradial.pfixmax);
% result.error=zeros(1,paramradial.nbfoncmax);
% result.errorres=zeros(1,paramradial.nbfoncmax);
% result.errorrayg=zeros(1,paramradial.nbfoncmax);
% result.errorraygi=zeros(1,paramradial.nbfoncmax);
% result.time = zeros(1,paramradial.nbfoncmax);
% result.rayg=cell(1,paramradial.nbfoncmax);
% result.Rayg=cell(1,paramradial.nbfoncmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;


%% D�finition des normes
Am = expect(A);
Am = 1/2*(Am+Am');
Am = speye(n);


if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

%% D�finition des solveurs
toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');

if directsolve
    localstosolver = @(aU,fU) solve(aU,fU);
else
    localstosolver = @(aU,fU) cgs(aU,fU,toliter,[],'noupdate');
end

clock0 = clock;


%% R�utilisation
if (paramradial.reuse && nargin>=4 && ~isempty(ureuse))
    [u,flagreuse,resultini,bini]=solve_reuse(GSD,A,b,ureuse,varargin{:});
else
    flagreuse = true;
end


if flagreuse
    u = PCRADIALMATRIX([n,1],PC);
    bini = b;
    result.raygini=0;
    result.errorini=1;
else
    result.raygini=resultini.rayg;
    result.Raygini=resultini.Rayg;
    
    switch errorindicator
        case 'rayleigh'
            result.errorini=1;
        case 'reference'
            result.errorini=mynorm(u-uref)/mynorm(uref);
        case 'residual'
            result.errorini=mynorm(bini)/normb;
    end
    if display_
        fprintf('  Initialisation : %3d fonctions',length(u))
    end
    if ~strcmp(errorindicator,'rayleigh')
        if display_
            fprintf(' - erreur  : %3d \n',result.errorini)
        end
        if result.errorini<getparam(GSD,'tolini')
            fprintf('  Convergence apres initialisation - erreur : %3d\n',result.errorini)
            return
        end
    end
end

result.nfonctions = length(u);

bu = bini;

result.totaltime=etime(clock,clock0);


%% Construction des couples de fonctions
nbfoncmax=paramradial.nbfoncmax;
for j=1:1
    
    
    inittype = getparam(GSD,'inittype');
    switch inittype
        case {'arnoldi','one'}
            
            result.nfonctions = result.nfonctions  + nbfoncmax  ;
            l0 = rand(PC);
            % l0 = ones(length(PC),1);l0=PCMATRIX(l0,[1,1],PC);
            V = zeros(n,0);
            if display_
                fprintf('initialisation arnoldi\n')
            end
            for ii=1:nbfoncmax
                f0=expectmtimes(bu,l0);
                L0 = expectmtimes(A,l0,l0);
                U = L0\f0;
                H0=prodscal(U,U,Am);
                for jjj=1:size(V,2)
                    U=U-V(:,jjj)*(V(:,jjj)'*Am*U);
                end
                H=sqrt(U'*Am*U);
                if display_;fprintf('fonction %d , residu = %.3e\n',ii,full(H/H0));end
                if H/H0<orthocrit
                    if display_
                        disp(' -> break')
                    end
                    break
                end
                U = U / H ;
                V=[V,U];
                fU = U'*bu;
                aU = U'*A*U;
                [l,flag] = localstosolver(aU,fU);
                l0 = l ; U0 = U ;
            end
            
            fU = V'*bu;
            aU = V'*A*V;
            [l,flag]=localstosolver(aU,fU);
            l = expand(l);
            V0=V;
            l0 = l;
            Rayg = double(full(expectmtimes(l,fU')));
            rayg = trace(Rayg);
            
            if ischarin('reference',varargin)
                u = PCRADIALMATRIX(V,[n,1],l);
                result.erroriter(1) =  norm(u-uref)/norm(uref);
                result.erroriter_A(1) =  norm(u-expand(uref),A)/norm(expand(uref),A);
                
                fprintf('  iteration %3d -> erreur   reference = %.3e  \r',0,result.erroriter(1));
                fprintf('  iteration %3d -> erreur-A reference = %.3e  \r',0,result.erroriter_A(1));
            end
            
            
        case 'random'
            if display_
                fprintf('initialisation random lambda\n')
            end
            l = rand(nbfoncmax,length(PC));l=PCMATRIX(l,[nbfoncmax,1],PC);
            Rayg = zeros(nbfoncmax);
            rayg = 0;
            
            
            if ischarin('reference',varargin) | paramradial.pfixmax==0
                m = numel(l);
                
                L0 = expect(A,l,l);
                if m>1
                    L0 = assembleblock(L0);
                end
                
                f0 = expect(bu,l);
                if m>1
                    f0 = assembleblock(f0);
                end
                V = L0\f0;
                V = reshape(V,[n,m]);
                u = PCRADIALMATRIX(V,[n,1],l);
                if ischarin('reference',varargin)
                    result.erroriter(1) =  norm(u-uref)/norm(uref);
                    result.erroriter_A(1) =  norm(u-expand(uref),A)/norm(expand(uref),A);
                    fprintf('  iteration %3d -> erreur   reference = %.3e  \r',0,result.erroriter(1));
                    fprintf('  iteration %3d -> erreur-A reference = %.3e  \r',0,result.erroriter_A(1));
                end
                
                V0 = V ;
            end
            
            
            
    end
    
    
    result.Rayg{j}{1}=Rayg;
    result.rayg{j}(1)=rayg;
    
    try
        result.stoV{1} = V0;
    end
    result.stoL{1} = l;
    
    
    for kk=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE
        
        m = numel(l);
        
        L0 = expect(A,l,l);
        if m>1
            L0 = assembleblock(L0);
        end
        
        f0 = expect(bu,l);
        if m>1
            f0 = assembleblock(f0);
        end
        V = L0\f0;
        V = reshape(V,[n,m]);
        
        V = orth(full(V));
        
        if size(V,2)<nbfoncmax
            fprintf('\nwarning : non linearly independent space functions  \n -> subspace dimension set to %d\n',size(V,2))
        end
        
        if ischarin('reference',varargin) & isa(uref,'PCRADIALMATRIX')
            result.subspaceangle_ref(kk)=max(subspaceangle(full(V),full(double(getV(uref)))));
            fprintf('sinus of the largest principal angle between iterate and reference = %d\n',...
                sin(max(result.subspaceangle_ref(kk))))
        end
        
        if ~(kk==1 & strcmp(inittype,'random'))
            result.subspaceangle_conv(kk)=max(subspaceangle(full(V),full(V0)));
            fprintf('sinus of the largest principal angle between iterates = %d\n',...
                sin(max(result.subspaceangle_conv(kk))))
        end
        
        
        result.stoV{kk+1} = V;
        
        V0=V;
        
        fU = V'*bu;
        aU = V'*A*V;
        [l,flag]=localstosolver(aU,fU);
        l = expand(l);
        
        result.stoL{kk+1}= l ;
        
        Rayg = double(full(expectmtimes(l,fU')));
        rayg = trace(Rayg);
        
        result.Rayg{j}{kk+1}=Rayg;
        result.rayg{j}(kk+1)=rayg;
        errorpf(kk)=abs((result.rayg{j}(kk+1)-result.rayg{j}(kk))/result.rayg{j}(kk+1));
        
        if display_
            fprintf('  iteration %3d -> erreur = %.3e  \r',kk,errorpf(kk));
        end
        
        
        if ischarin('reference',varargin)
            u = PCRADIALMATRIX(V,[n,1],l);
            result.erroriter(end+1) =  norm(u-uref)/norm(uref);
            result.erroriter_A(end+1) =  norm(u-expand(uref),A)/norm(expand(uref));
            
            fprintf('  iteration %3d -> erreur   reference = %.3e  \r',kk,result.erroriter(end));
            fprintf('  iteration %3d -> erreur-A reference = %.3e  \r',kk,result.erroriter_A(end));
            
        end
        
        
        if errorpf(kk)<paramradial.pfixtol
            break
        end
        
    end
    
    
    if j==1
        result.errorrayg(j) = sqrt(abs(result.rayg{j}-result.raygini)/((result.rayg{j})));
    else
        result.errorrayg(j)=sqrt(abs(result.rayg{j}-result.rayg{j-1})/(result.rayg{j}));
    end
    
    
    u = PCRADIALMATRIX(V,[n,1],l);
    bu = bu - A*u;
    
    % CALCUL ERREUR
    switch errorindicator
        case 'residual'
            result.error(j) = mynorm(bu)/normb;
        case 'rayleigh'
            result.error(j)= result.errorrayg(j)  ;
        case 'reference'
            result.error(j)=mynorm(u-uref)/mynorm(uref);
    end
    if ischarin('testerror',varargin)
        % result.errorA(j)=mynorm(u-uref,A)/mynorm(uref,A);
        result.errorres(j)=mynorm(bu)/normb;
        
    end
    
    if display_
        fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        if ischarin('testerror',varargin)
            fprintf('  %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
        end
    end
    
    result.addedfunctions = numel(l);
    result.time(j)=etime(clock,clock0);
    
    if result.error(j)<paramradial.tol
        hasconv=true;
        fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        break
    end
    
end


hasconv=false;

result.totaltime=etime(clock,clock0);


%% Finalisation : decomposition spectrale de la solution

if getparam(GSD,'finalunique')
    tolfactSD = getparam(GSD,'finaluniquefacttol');
    u=uniqueV(u,paramradial.tol*tolfactSD);
    
elseif getparam(GSD,'finalSD')
    tolfactSD = getparam(GSD,'finalSDfacttol');
    u = spectral_decomposition(u,'tol',paramradial.tol*tolfactSD,...
        'nbfoncmax',getm(u));
    
end

if getparam(GSD,'finalupdate') %%%% FINAL UPDATING
    U = double(getV(u));
    fU = U'*b;
    aU = U'*A*U;
    [l,flag]=localstosolver(aU,fU);
    
    u = PCRADIALMATRIX(U,[n,1],l);
    bu = b - A*u;
    switch errorindicator
        case 'reference'
            result.error(j)=mynorm(u-uref)/mynorm(uref);
        otherwise
            result.error(j) = mynorm(bu)/normb;
    end
    
    disp('  Reactualisation')
    fprintf('  nb fonctions %3d - erreur  : %3d \n',result.nfonctions,result.error(j))
end %%%% END FINAL UPDATING

%% Verdict sur la convergence
if ~hasconv
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
end

if display_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end

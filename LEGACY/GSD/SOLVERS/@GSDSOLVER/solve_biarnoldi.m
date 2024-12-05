function [u,result] = solve_arnoldi(GSD,A,b,ureuse,varargin)
% function u = solve_arnoldi(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% algorithme de type arnoldi
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

display_ = getparam(GSD,'display');
fprintf('GSD solveur arnoldi ... ')
if display_
    fprintf('\n')
end
uref = getcharin('reference',varargin);
paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.nbfoncmaxsimul = getparam(GSD,'nbfoncmaxsimul');
orthocrit = getparam(GSD,'orthocrit');
restart = getparam(GSD,'restart');
paramradial.tol = getparam(GSD,'tol');
paramradial.reuse = getparam(GSD,'reuse');
errorindicator = getparam(GSD,'errorindicator');
update = getparam(GSD,'update');
% deflation = getparam(GSD,'deflation');


if ~isempty(uref)
    errorindicator='reference';
end


[A,b,PC]=pcsystemupdate(A,b,varargin{:});
n = size(A,1);
result.nfonctions = 0 ;
result.addedfunctions = 0;



clock0 = clock;

%% definition des normes et solveurs

Am = speye(n);

if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');

if directsolve
    localstosolver = @(aU,fU) solve(aU,fU);
else
    localstosolver = @(aU,fU) cgs(aU,fU,toliter,[],'noupdate');
end


switch errorindicator
    case 'rayleigh'
        result.errorini=1;
    case 'reference'
        result.errorini=mynorm(u-uref)/mynorm(uref);
    case 'residual'
        result.errorini=mynorm(bini)/normb;
    case 'none'
        result.errorini=1;
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

result.nfonctions = getm(u);

uini = u ;
bu = bini;

result.addedfunctions = 0;


Vu = sparse(n,0);

l0 = zeros(1,PC);
l = zeros(1,PC);
U = zeros(n,1);
U0 = zeros(n,1);
stoU = sparse(n,0);
stol = zeros(0,1,PC);

clock0=clock;

%% Construction des couples de fonctions
for k=0:restart %%% RESTARTS DE ARNOLDI
    if result.addedfunctions>=paramradial.nbfoncmax
        break
    end
    if k>=1
        fprintf('Actualisation #%d \n',k)
    end
    
    if update==0 || k==0
        V = zeros(n,0);
        L = zeros(0,1,PC);
    end
    
    switch getparam(GSD,'inittype')
        case 'one'
            l0=one(PC);
        case 'random'
            l0 = rand(PC);
        case 'allone'
            l0 = ones(length(PC),1);
            l0 = PCMATRIX(l0,[1,1],PC);
            
    end
    
    alpha0=norm(l0);l0=l0/alpha0;
    
    
    for i=1:paramradial.nbfoncmaxsimul %%% CONSTRUCTION DE LA BASE DE KRYLOV
        f0=expectmtimes(bu,l0);
        L0 = expectmtimes(A,l0,l0);
        
        
        U = L0\f0;
        
        H0=prodscal(U,U,Am);
        
        for j=1:size(V,2)
            U=U-V(:,j)*(V(:,j)'*Am*U);
        end
        H=sqrt(U'*Am*U);
        
        if display_
            fprintf('residu = %.3e',full(H/H0))
        end
        
        if H/H0<orthocrit
            if display_
                disp(' -> break')
            end
            break
        end
        
        
        U = U / H ;
        
        V=[V,U];
        
        result.addedfunctions=result.addedfunctions+1;
        if result.addedfunctions>=paramradial.nbfoncmax
            fprintf('\n')
            break
        end
        
        fprintf(' -> %d fonctions\n',size(V,2))
        
        fU = U'*bu;
        aU = U'*A*U;
        
        [l,flag] = localstosolver(aU,fU);
        
        l0 = l ; U0 = U ;
        
        result.time(result.addedfunctions)=etime(clock,clock0);
        switch errorindicator
            case 'reference'
                fUtemp = V'*bu;
                aUtemp = V'*A*V;
                [ltemp,flag] = localstosolver(aUtemp,fUtemp);
                uadd = PCRADIALMATRIX(V,[n,1],ltemp);
                utemp = u + uadd;
                ii=result.addedfunctions;
                result.error(ii)=mynorm(utemp-uref)/mynorm(uref);
        end
    end
    
    
    if update
        buup = bini;
    else
        buup = bu;
    end
    
    fU = V'*buup;
    aU = V'*A*V;
    [l,flag] = localstosolver(aU,fU);
    
    result.Rayg{i,k+1} = expectmtimes(l,fU');
    result.rayg{i,k+1} = trace(result.Rayg{i,k+1});
    
    if getparam(GSD,'subspaceiteration')>0
        bsub = buup;
        for kk=1:getparam(GSD,'subspaceiteration')
            if display_
                fprintf('subspace iteration %d \n',kk)
            end
            
            subspaceiteration
            
            
            result.raygsub{k+1}(kk)=rayg;
            errorpf(kk)=abs((result.raygsub{k+1}(kk)-(kk>1)*result.raygsub{k+1}(kk-(kk>1)))/result.raygsub{k+1}(kk));
            
        end
        
    end
    
    
    if update
        uadd = PCRADIALMATRIX(V,[n,1],l);
        u = uini + PCRADIALMATRIX(V,[n,1],l);
        bu = bini - A*uadd;
    else
        uadd = PCRADIALMATRIX(V,[n,1],l);
        u = u + uadd;
        bu = bu - A*uadd;
    end
    
    if getparam(GSD,'righthandSD');
        bu = righthandSD(GSD,expand(bu));
    end
    
    
    ii=result.addedfunctions;
    
    result.time(ii)=etime(clock,clock0);
    
    switch errorindicator
        case 'residual'
            result.error(ii) = mynorm(bu)/normb;
        case 'reference'
            result.error(ii)=mynorm(u-uref)/mynorm(uref);
        case 'none'
            result.error(ii)=1;
        otherwise
            result.error(ii)=1;
            warning('type d''erreur pas bon')
    end
    result.nfonctions = length(u);
    
    result.time(ii)=etime(clock,clock0);
    
    if result.error(ii)<paramradial.tol
        hasconv=true;
        fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
        break
    end
    hasconv=false;
    if display_
        fprintf('%d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
    end
    
end

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


%% Verdict sur la convergence
if ~hasconv && paramradial.nbfoncmax>0
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
end

if display_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end



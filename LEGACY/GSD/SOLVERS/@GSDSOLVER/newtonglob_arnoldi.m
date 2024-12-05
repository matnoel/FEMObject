function [u,result] = newtonglob_arnoldi(GSD,N,S,b,X,varargin)


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

if ~isempty(uref)
    errorindicator='reference';
end

PC = getPC(X);
P = getP(PC);
n=size(b,1);

result.nfonctions = 0 ;
result.addedfunctions = 0;


clock0 = clock;

%% definition des normes et solveurs

if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');


result.addedfunctions = 0;

Vu = sparse(n,0);

l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
stoU = sparse(n,0);
U0 = sparse(n,1);

clock0=clock;


N = setparam(N,'increment',true)
N = setparam(N,'stopini',false);

u = PCRADIALMATRIX([n,1],PC);

%% Construction des couples de fonctions
for k=0:restart %%% RESTARTS DE ARNOLDI
    if result.addedfunctions>=paramradial.nbfoncmax
        break
    end
    if k>=1
        fprintf('Actualisation #%d \n',k)
    end
    
    
    V = sparse(n,0);
    l0 = initlambda(GSD,PC);
    
    alpha0=norm(l0);l0=l0/alpha0;
    
    
    for i=1:paramradial.nbfoncmaxsimul %%% CONSTRUCTION DE LA BASE DE KRYLOV
        f0=expectmtimes(b,l0);
        
        U = solve(N,f0,@(U) expectmtimes(calc_fintpc(S,(u+U*l0),X),l0),...
            @(U) expectmtimes(calc_rigitangpc(S,(u+U*l0),X),l0,l0),zeros(n,1));
        
        H0=prodscal(U,U);
        
        for j=1:size(V,2)
            U=U-V(:,j)*(V(:,j)'*U);
        end
        
        H=sqrt(U'*U);
        
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
        
        fU = U'*b;
        
        l = solve(N,fU,@(l) U'*calc_fintpc(S,(u+U*l),X),...
            @(l) U'*calc_rigitangpc(S,(u+U*l),X)*U,zeros(1,PC));
        
        %lr = rand(PC);lr=lr/norm(lr)*norm(l);
        %l = l + lr*0.2;
        l0 = l ; U0 = U ;
        
        result.time(result.addedfunctions)=etime(clock,clock0);
        
    end
    
    fU = V'*b;
    l0 = zeros(size(V,2),X);
    l = solve(setparam(N,'tol',paramradial.tol/10),...
        fU,...
        @(l) V'*calc_fintpc(S,u+PCRADIALMATRIX(V,[n,1],expand(l)),X),...
        @(l) V'*calc_rigitangpc(S,u+PCRADIALMATRIX(V,[n,1],expand(l)),X)*V,l0);
    
    result.Rayg{i,k+1} = expectmtimes(l,fU');
    result.rayg{i,k+1} = trace(result.Rayg{i,k+1});
    
    uadd = PCRADIALMATRIX(V,[n,1],l);
    u=u+uadd ;
    
    bu = b - calc_fintpc(S,u,X) ;
    
    Vu = [Vu,V];
    
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
            error('type d''erreur pas bon')
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
if ~hasconv
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
end

if display_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end



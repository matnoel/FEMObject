function [u,result] = solve_arnoldi_separated(GSD,A,b,ureuse,varargin)
% function u = solve_arnoldi_separated(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

display_ = getparam(GSD,'display');

fprintf('\n GSD solveur ... ')
if display_
    fprintf('\n')
end

paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
update = getparam(GSD,'update');
update=true;
paramradial.nbfoncmaxsimul = getparam(GSD,'nbfoncmaxsimul');
orthocrit = getparam(GSD,'orthocrit');
restart = getparam(GSD,'restart');
paramradial.orthoduringpfix = getparam(GSD,'orthoduringpfix');

errorindicator = getparam(GSD,'errorindicator');
uref = getcharin('reference',varargin);

if ~isempty(uref)
    errorindicator='reference';
end

n=size(A,1);

errorpf=zeros(1,paramradial.pfixmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;

PC = getPC(A);


%% D�finition des normes
Am = speye(n);
if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

%% D�finition des solveurs

GSD = updatelocalstosolver(GSD);
localstosolveriter = getparam(GSD,'localstosolveriter');
localstosolverupdate = getparam(GSD,'localstosolverupdate');
toliter = getparam(GSD,'toliter');
tolupdate = getparam(GSD,'tolupdate');

% localstosolver = getparam(GSD,'localstosolver');

localstosolveriter = @(A,b) solve(localstosolveriter,A,b,toliter);
localstosolverupdate = @(A,b)  solve(localstosolverupdate,A,b,tolupdate);

clock0 = clock;

result.totaltime=etime(clock,clock0);

u = PCTPRADIALMATRIX([n,1],PC);
result.raygini=0;
result.errorini=1;
result.Raygini = 1;

l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
stoU = sparse(n,0);
stol = cell(1,0);
U0 = sparse(n,1);

nbfoncmax = paramradial.nbfoncmax;
result.nfonctions = 0;
result.addedfunctions = 0;


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
    end
    
    switch getparam(GSD,'inittype')
        case 'one'
            l0 = one(PC);
        case 'random'
            l0 = rand(PC);
        otherwise
            error('pas defini')
            
    end
    
    alpha0=norm(l0);l0=l0/alpha0;
    
    for i=1:paramradial.nbfoncmaxsimul
        
        f0=expect(b,l0);
        for kk=1:getm(u)
            f0 = f0 - expect(A*getV(u,kk),l0,getL(u,kk));
        end
        L0 = expect(A,l0,l0);
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
        
        
        fU = U'*b;
        for kk=1:getm(u)
            fU = fU - (U'*A*getV(u,kk))*getL(u,kk);
        end
        
        aU = U'*A*U;
        
        [l,flag]=localstosolveriter(aU,fU);
        l0 = l ; U0 = U ;
        result.time(result.addedfunctions)=etime(clock,clock0);
        
    end
    
    
    fU = V'*b;
    aU = V'*A*V;
    
    if display_
        fprintf('Actualisation des variables aleatoires\n')
    end
    
    [l,flag]=localstosolverupdate(aU,fU);
    
    result.Rayg{i,k+1} = double(full(expectmtimes(l,fU')));
    result.rayg{i,k+1} = trace(result.Rayg{i,k+1});
    
    u = PCTPRADIALMATRIX(V,[n,1],l);
    
    ii=result.addedfunctions;
    
    result.time(ii)=etime(clock,clock0);
    
    
    switch errorindicator
        case 'residual'
            result.error(ii) = mynorm(b-A*u)/normb;
        case 'reference'
            result.error(ii)=mynorm(u-uref)/mynorm(uref);
        case 'none'
            result.error(ii)=1;
        otherwise
            result.error(ii)=1;
            warning('type d''erreur pas bon')
    end
    
    result.nfonctions = getm(u);
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
    
    
end %%%%%%%%%%%%% FIN RESTART


result.totaltime=etime(clock,clock0);



%% Verdict sur la convergence
if ~hasconv && paramradial.nbfoncmax>0
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
end

if display_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end
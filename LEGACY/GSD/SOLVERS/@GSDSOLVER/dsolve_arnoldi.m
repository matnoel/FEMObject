function [u,result] = dsolve_arnoldi(GSD,T,b,A,B,u0,ureuse,varargin)
% function [u,result] = dsolve_arnoldi(GSD,T,b,A,B,u0,ureuse,varargin)


if nargin>=6 && ~isempty(u0)
    warning('verifier prise en compte de la condition initiale')
    b0 = create_initial_condition(T,A*u0);
    b = b+b0;
end

fprintf('\nGSD solveur arnoldi ... \n')

paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.nbfoncmaxsimul = getparam(GSD,'nbfoncmaxsimul');
orthocrit = getparam(GSD,'orthocrit');
restart = getparam(GSD,'restart');
paramradial.tol = getparam(GSD,'tol');
paramradial.reuse = getparam(GSD,'reuse');
paramradial.orthogonalizeLresidual = getparam(GSD,'orthogonalizeLresidual');
errorindicator = getparam(GSD,'errorindicator');
update = getparam(GSD,'update');


uref = getcharin('reference',varargin);
varargin = delcharin('reference',varargin);
if ~isempty(uref)
    errorindicator='reference';
end


toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');
dislay_ = getparam(GSD,'display');


[b,A,B] = init_resolution(T,b,A,B);
[b,A,B,PC]=pcsystemupdate(b,A,B,varargin{:});

n = size(A,1);
nt = length(T);
P=length(PC);

%% definition des normes et solveurs

if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end



%% R�utilisation
if (paramradial.reuse && nargin>=7 && ~isempty(ureuse))
    [u,flagreuse,resultini,bini]=dsolve_reuse(GSD,T,A,B,b,ureuse,varargin{:});
else
    flagreuse = true;
end


if flagreuse
    u = PCTIMEMATRIX(PCRADIALMATRIX([n,nt],PC),T,[n,1]);
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
        case 'none'
            result.errorini=1;
    end
    if dislay_
        fprintf('  Initialisation : %3d fonctions',getm(u))
    end
    if ~strcmp(errorindicator,'rayleigh')
        if dislay_
            fprintf(' - erreur  : %3d \n',result.errorini)
        end
        if result.errorini<getparam(GSD,'tolini')
            fprintf('  Convergence apres initialisation - erreur : %3d\n',result.errorini)
            return
        end
    end
end

bu = bini;
uini = u;

result.nfonctions = getm(u);
result.addedfunctions = 0;


clock0=clock;


for k=0:restart
    if result.addedfunctions>=paramradial.nbfoncmax
        break
    end
    if k>=1
        fprintf('Actualisation #%d \n',k)
    end
    
    if update==0 || k==0
        V = cell(1,0);
        DV = cell(1,0);
    end
    l0 = rand(PC);
    alpha0=norm(l0);l0=l0/alpha0;
    
    
    for i=1:paramradial.nbfoncmaxsimul
        f0=expectmtimes(bu,l0);
        A0 = expectmtimes(A,l0,l0);
        B0 = expectmtimes(B,l0,l0);
        U = dsolve(T,f0,A0,B0);
        
        H0=norm(U);
        
        for j=1:size(V,2)
            U=U-V{j}*prodscal(V{j},U);
        end
        
        H=norm(U);
        
        if dislay_
            fprintf('residu = %.3e',full(H/H0))
        end
        
        if H/H0<orthocrit
            if dislay_
                disp(' -> break')
            end
            break
        end
        
        
        U = U / H ;
        
        DU = diff(U);
        
        V=[V,{U}];
        DV=[DV,{DU}];
        
        
        result.addedfunctions=result.addedfunctions+1;
        if result.addedfunctions>=paramradial.nbfoncmax
            fprintf('\n')
            break
        end
        
        if dislay_
            fprintf(' -> %d fonctions\n',size(V,2))
        end
        
        
        fU = integratemtimes(U',bu);
        aU = integratemtimes(U',A,DU);
        bU = integratemtimes(U',B,U);
        aU = aU+bU;
        result.time(result.addedfunctions)=etime(clock,clock0);
        
        if directsolve
            l = solve(aU,fU);
        else
            [l,flag] = cgs(aU,fU,toliter,[],'noupdate');
        end
        
        l0 = l ; U0 = U ;
        
    end
    clear aU bU fU;
    
    if update
        buup = bini;
    else
        buup = bu;
    end
    
    for i=1:length(V)
        fU{i} = integratemtimes(V{i}',buup);
        for j=1:length(V)
            aU{i}{j} = integratemtimes(V{i}',A,DV{j});
            bU{i}{j} = integratemtimes(V{i}',B,V{j});
            aU{i}{j} = aU{i}{j}+bU{i}{j};
        end
        aU{i} = horzcat(aU{i}{:});
    end
    fU = vertcat(fU{:});
    aU = vertcat(aU{:});
    
    if directsolve
        l = solve(expand(aU),fU);
    else
        [l,flag] = cgs(aU,fU,toliter,[],'noupdate');
    end
    
    result.Rayg{i,k+1} = double(cell2mat(full(expect(fU,l))));
    result.rayg{i,k+1} = trace(result.Rayg{i,k+1});
    
    for i=1:length(V)
        Vi = PCRADIALMATRIX(getvalue(V{i}),[n,nt],l(i));
        if update && i==1
            u = uini + PCTIMEMATRIX(Vi,T,[n,1]);
        else
            u = u + PCTIMEMATRIX(Vi,T,[n,1]);
        end
        buup = buup - (A*DV{i}+B*V{i})*l(i);
        
    end
    bu=buup;
    
    if paramradial.orthogonalizeLresidual && ispcradial(bu)
        if dislay_
            fprintf('Sorting radial functions of residual : %d -> ',getm(bu))
        end
        bu = orthogonalizeL(bu,paramradial.tol);
        if dislay_
            fprintf('%d functions\n',getm(bu))
        end
    end
    
    
    result.nfonctions = length(getvalue(u));
    
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
    
    result.time(ii)=etime(clock,clock0);
    
    if result.error(ii)<paramradial.tol
        hasconv=true;
        fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
        break
    end
    hasconv=false;
    if dislay_
        fprintf('%d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
    end
    
end

result.totaltime=etime(clock,clock0);


if ~hasconv
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(ii))
end

fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))


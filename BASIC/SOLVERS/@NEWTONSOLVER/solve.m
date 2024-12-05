function [u,result] = solve(N,b,calcAu,calcAtang,u0,Atang,solver,S)
% function u = solve(N,b,calcAu,calcAtang,u0,Atang,solver,S)
% solveur Newton : resolution de A(u)=b
% N : NEWTONSOLVER
% b : second membre
% calcAu : fonction telle que calcAu(u) = A(u)
% calcAtang : fonction telle que calcAtang(u) = A'(u) (matrice tangente en u)
% u0 : initialisation
% Atang : matrice forcee dans Newton constant
% solver : fonction donnant le solveur des systemes lineaires
%          solver(A,b) resout Au=b
% S : (optionel) maillage associe a la solution u

tic;

type = getparam(N,'type');
display_ = getparam(N,'display');
if display_
    fprintf('\n ---------------------------------------');
    fprintf('\n ------------ Newton solver ------------')
    fprintf('\n ---------------------------------------\n');
end

calcAu = fcnchk(calcAu);

try
    calcAtang = fcnchk(calcAtang);
catch
    warning('pas de fonction donnant la matrice tangente')
end

tol = getparam(N,'tol');
tolstagn = getparam(N,'tolstagn');
tolreact = getparam(N,'tolreact');
maxiter = getparam(N,'maxiter');
increment = getparam(N,'increment');

if nargin<5 || isempty(u0)
    u0 = zeros(size(b));
end

if nargin<6 || isempty(Atang) || ~(strcmpi(type,'constant') || strcmpi(type,'manual'))
    Atang = calcAtang(u0);
end

if nargin==7
    solver = fcnchk(solver);
else
    solver = @solve;
end

r0 = b-calcAu(u0);

if norm(b)<eps
    nb = norm(calcAu(u0));
else
    nb = norm(b);
end

err = norm(r0)/nb;
result.errorini = err;

if err<tol
    if getparam(N,'stopini')
        if display_
            fprintf('\nNewton converged at initialization with residual error %.2d\n',err);
        end
        u = u0;
        return
    else
        if display_
            fprintf('\nNewton converged at initialization with residual error %.2d -> NO STOP \n',err);
        end
    end
end

if display_
    fprintf('\nNewton initialization : residual = %.2d',err);
end

for i=1:maxiter
    
    if increment
        if nargin(solver)==3
            [du,result.solver{i}] = solver(Atang,r0,u0);
        else
            [du,result.solver{i}] = solver(Atang,r0);
        end
        
        if nargin==8
            u = unfreevector(S,u0)+unfreevector(S,du);
        else
            u = du+u0;
        end
    else
        if nargin==8
            warning('ne marche pas avec dirichlet non-homogene')
        end
        if nargin(solver)==3
            [u,result.solver{i}] = solver(Atang,r0+Atang*u0,u0);
        else
            [u,result.solver{i}] = solver(Atang,r0+Atang*u0);
        end
    end
    %%%% le mettre dans le solver GSDSOLVER comme option et le faire
    % a la fin des solve_power, solve_arnoldi ...
    % if isa(u,'PCRADIALMATRIX')
    %     u = spectral_decomposition(u,'tol',min(tol,1e-10),'nbfoncmax',length(u));
    %     u = uniqueV(u);
    % end
    
    r = b-calcAu(u);
    err = norm(r)/nb;
    result.error(i) = err;
    result.time(i) = toc;
    if err<tol
        break
    else
        if display_
            fprintf('\nNewton iteration #%3.d : residual = %.2d',i,err);
        end
    end
    
    errstagn = norm(u0-u)/norm(u+u0);
    
    if (strcmpi(type,'tangent') || strcmpi(type,'full') || errstagn<tolstagn || err>tolreact) && ~(strcmpi(type,'constant') || strcmpi(type,'manual'))
        if display_
            fprintf(' (update tangent)')
        end
        Atang = calcAtang(u);
    end
    
    u0 = u;
    r0 = r;
end


if err<tol
    if display_
        fprintf('\nNewton converged at iteration #%d with residual error %.2d\n',i,err);
    end
else
    fprintf('\nNewton stopped at iteration #%d with residual error %.2d\n',i,err);
end

if display_
    fprintf('Elapsed time = %f s\n',toc)
end

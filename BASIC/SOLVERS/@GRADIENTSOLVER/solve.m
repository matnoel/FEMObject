function [u,result] = solve(G,b,calcAu,calcBu,u0)
% function u = solve(G,b,calcAu,calcBu,u0,solveur)
% solveur : resolution de A(u)+B(u)=b
% avec B lineaire et phi nonlineaire
% G : GRADIENT
% b : second membre
% calcAu : fonction telle que calcAu(u) = A(u)
% calcBu : fonction telle que calcBu(u) = B(u)
% u0 : initialisation

tic;

display_ = getparam(G,'display');
if display_
    fprintf('\n -----------------------------------------------------');
    fprintf('\n ------------ Alternated direction solver ------------')
    fprintf('\n -----------------------------------------------------\n');
end

calcAu = fcnchk(calcAu);
calcBu = fcnchk(calcBu);

tol = getparam(G,'tol');
maxiter = getparam(G,'maxiter');

if nargin<5 || isempty(u0)
    u0 = zeros(size(b));
end

r0 = b-calcAu(u0)-calcBu(u0);

if norm(b)<eps
    nb = norm(calcAu(u0)+calcBu(u0));
else
    nb = norm(b);
end

err = norm(r0)/nb;
result.errorini = err;

if err<tol
    if getparam(G,'stopini')
        if display_
            fprintf('\nAlt dir converged at initialization with residual error %.2d\n',err);
        end
        u = u0;
        return
    else
        if display_
            fprintf('\nAlt dir converged at initialization with residual error %.2d -> NO STOP \n',err);
        end
    end
end

if display_
    fprintf('\nAlt dir initialization : residual = %d',err);
end

errstagn = 1;
for i=1:maxiter
    r = b-calcAu(u);
    err = norm(r)/nb;
    result.error(i) = err;
    result.time(i) = toc;
    if err<tol
        break
    else
        if display_
            fprintf('\nAlt dir iteration #%3.d : residual = %.2d',i,err);
        end
    end
    
    errstagn = norm(u0-u)/norm(u+u0);
    
    u0 = u;
    r0 = r;
end


if err<tol
    if display_
        fprintf('\nAlt dir converged at iteration #%d with residual error %.2d\n',i,err);
    end
else
    fprintf('\nAlt dir stopped at iteration #%d with residual error %.2d\n',i,err);
end

if display_
    fprintf('Elapsed time %f s\n',toc)
end

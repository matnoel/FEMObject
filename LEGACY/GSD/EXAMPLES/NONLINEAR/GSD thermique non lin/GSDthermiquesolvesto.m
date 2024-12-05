function l=GSDthermiquesolvesto(PC,tol,a1,a2,a3,b,varargin)


display_ = ischarin('display',varargin);
if display_
    fprintf('\n RESOLUTION PROB STOCHASTIQUE\n')
end

l0 = zeros(PC,1);
l=l0;
for k=1:100
    res = b - a1*l - a2*l*l - a3*l*l*l;
    nres = norm(res)/norm(b);
    if nres<tol
        break
    end
    if display_
        fprintf('iteration %d , residual = %.2d\n',k,nres);
    end

    AT = a1+2*a2*l + 3*a3*l*l;
    dl = solve(AT,res);
    l = l+dl;

    l0=l;
end

if ~display_
    fprintf('(fM(U) %d iterations)',k);
end

if nres>=tol
    warning('NON CONVERGENCE');
end

function at = dderivative(L,ut,vt,b,M,K,C,solver)
% at = dderivative(L,ut,vt,b,M,K,C,solver)
% Computes the second order derivative a=u'' of u 
% by solving the differential equation Mu''+Cu'+Ku=b
% L: DGTIMESOLVER
% b: right-hand side
% M,K,C: double or random matrices
% solver: function for solving linear systems
%         solver(A,b) solves Au=b
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'
% at: TIMEMATRIX of acceleration a=v'=u''

n = size(M,1);
if nargin<7 || isempty(C)
    C = sparse(n,n);
end

if nargin>=8 && ~isempty(solver)
    solver =  fcnchk(solver);
else
    if ~istime(M) && isa(M,'double') && getparam(L,'lu')
        [LM,UM] = lu(M);
        solver = @(M,b) solve(UM,(LM\b));
    else
        solver = @solve;
    end
end

T = gettimemodel(L);
nt = getnt(T);

at = cell(1,nt);
for i=1:nt
    i0 = 2*i-1; % t0+
    i1 = 2*i; % t1-
    b0 = b{i0};
    b1 = b{i1};
    if istime(M)
        M0 = getmatrixatstep(M,i0);
        M1 = getmatrixatstep(M,i1);
    else
        M0 = M;
        M1 = M;
    end
    if istime(C)
        C0 = getmatrixatstep(C,i0);
        C1 = getmatrixatstep(C,i1);
    else
        C0 = C;
        C1 = C;
    end
    if istime(K)
        K0 = getmatrixatstep(K,i0);
        K1 = getmatrixatstep(K,i1);
    else
        K0 = K;
        K1 = K;
    end
    at{i0} = solver(M0,b0-C0*vt{i0}-K0*ut{i0});
    at{i1} = solver(M1,b1-C1*vt{i1}-K1*ut{i1});
end

try
    at = horzcat(at{:});
catch
    warning('cell array')
end
at = TIMEMATRIX(at,T,[n,1]);

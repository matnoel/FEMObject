function vt = derivative(L,ut,b,A,B,solver)
% vt = derivative(L,ut,b,A,B,solver)
% Computes the first order derivative v of u 
% by solving the differential equation Au'+Bu=b
% L: DGTIMESOLVER
% b: right-hand side
% A,B: double or random matrices
% solver: function for solving linear systems
%         solver(A,b) solves Au=b
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'

n = size(A,1);
T = gettimemodel(L);
nt =  getnt(T);

if nargin>=6 && ~isempty(solver)
    solver =  fcnchk(solver);
else
    if ~istime(A) && isa(A,'double') && getparam(L,'lu')
        [LA,UA] = lu(A);
        solver = @(A,b) solve(UA,(LA\b));
    else
        solver = @solve;
    end
end

vt = cell(1,nt);
for i=1:nt
    i0 = 2*i-1; % t0+
    i1 = 2*i; % t1-
    b0 = b{i0};
    b1 = b{i1};
    if istime(A)
        A0 = getmatrixatstep(A,i0);
        A1 = getmatrixatstep(A,i1);
    else
        A0 = A;
        A1 = A;
    end
    if istime(B)
        B0 = getmatrixatstep(B,i0);
        B1 = getmatrixatstep(B,i1);
    else
        B0 = B;
        B1 = B;
    end
    vt{i0} = solver(A0,b0-B0*ut{i0});
    vt{i1} = solver(A1,b1-B1*ut{i1});
end

try
    vt = horzcat(vt{:});
catch
    warning('cell array')
end
vt = TIMEMATRIX(vt,T,[n,1]);

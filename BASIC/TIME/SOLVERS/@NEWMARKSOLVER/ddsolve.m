function [ut,result,vt,at] = ddsolve(L,b,M,K,C,u0,v0,solver)
% function [ut,result,vt,at] = ddsolve(L,b,M,K,C,u0,v0,solveur)
% Newmark solver: solves Mu''+Cu'+Ku=b
% L: NEWMARKSOLVER
% b: right-hand side
% M,K,C: double or random matrices
% solver: function for solving linear systems
%         solver(A,b) solves Au=b
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'
% at: TIMEMATRIX of acceleration a=v'=u''
% result: struct containing outputs
% result.totaltime : total CPU time

tStart = tic;

display_ = getparam(L,'display');
if display_
    fprintf('\n --------------------------------------------');
    fprintf('\n ------------ Newmark time solver -----------')
    fprintf('\n --------------------------------------------\n');
end

T = gettimemodel(L);
% t = gett(T);
nt = getnt(T);
dt = getdt(T);

alpha = getparam(L,'alpha');
beta = getparam(L,'beta');
gamma = getparam(L,'gamma');

if nargin<5
    C = [];
end

[b,M,K,C] = init_resolution(T,b,M,K,C);

if isuniform(T)
    A = M + (1-alpha)*beta*K*dt(1)^2;
    if ~isempty(C)
        A = A + (1-alpha)*gamma*C*dt(1);
    end
end

if nargin>=8 && ~isempty(solver)
    solver =  fcnchk(solver);
else
    if isuniform(T) && isa(A,'double')
        [LA,UA] = lu(A);
        solver = @(A,b) solve(UA,(LA\b));
    else
        solver = @solve;
    end
end

n = size(M,1);
ut = cell(1,length(T));
vt = cell(1,length(T));
at = cell(1,length(T));
if nargin<6 || isempty(u0)
    ut{1} = zeros(n,1);
else
    ut{1} = u0;
end
if nargin<7 || isempty(v0)
    vt{1} = zeros(n,1);
else
    vt{1} = v0;
end
b1 = b{1}-K*ut{1};
if ~isempty(C)
    b1 = b1-C*vt{1};
end
at{1} = solver(M,b1);

if display_
    fprintf('Newmark alpha=%.2f, beta=%.2f, gamma=%2f : ',alpha,beta,gamma);
end
for i=1:nt
    if display_
        pourcentage(i,nt)
    end
    
    ut_pred = ut{i} + vt{i}*dt(i) + (1/2-beta)*at{i}*dt(i)^2;
    vt_pred = vt{i} + (1-gamma)*at{i}*dt(i);
    bi = (1-alpha)*(b{i+1}-K*ut_pred) + alpha*(b{i}-K*ut{i});
    if ~isempty(C)
        bi = bi - (1-alpha)*C*vt_pred - alpha*C*vt{i};
    end
    if ~isuniform(T)
        A = M+(1-alpha)*beta*K*dt(i)^2;
        if ~isempty(C)
            A = A + (1-alpha)*gamma*C*dt(i);
        end
    end
    [at{i+1},result.solver{i}] = solver(A,bi);
    vt{i+1} = vt_pred + gamma*at{i+1}*dt(i);
    ut{i+1} = ut_pred + beta*at{i+1}*dt(i)^2;
    % vt{i+1} = vt{i} + ((1-gamma)*at{i} + gamma*at{i+1})*dt(i);
    % ut{i+1} = ut{i} + vt{i}*dt(i) + ((1/2-beta)*at{i} + beta*at{i+1})*dt(i)^2;
end

result.totaltime = toc(tStart);
if display_
    fprintf('Elapsed time = %f\n',result.totaltime)
end

try
    ut = horzcat(ut{:});
    vt = horzcat(vt{:});
    at = horzcat(at{:});
catch
    warning('cell array')
end

ut = TIMEMATRIX(ut,T,[n,1]);
vt = TIMEMATRIX(vt,T,[n,1]);
at = TIMEMATRIX(at,T,[n,1]);

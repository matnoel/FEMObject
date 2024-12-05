function [ut,result,vt] = dsolve(L,b,A,B,u0,solver)
% function [ut,result,vt] = dsolve(L,b,A,B,u0,solver)
% Euler solver: solves Au'+Bu=b
% L: EULERTIMESOLVER
% b: right-hand side
% A,B: double or random matrices
% u0: vector of initial conditions
% solver: function for solving linear systems
%         solver(A,b) solves Au=b
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'
% result: struct containing outputs
% result.totaltime : total CPU time

tStart = tic;

display_ = getparam(L,'display');
eulertype = getparam(L,'eulertype');
if display_
    fprintf('\n -------------------------------------------');
    fprintf('\n ------------ Euler time solver ------------')
    fprintf('\n -------------------------------------------\n');
end

T = gettimemodel(L);
% t = gett(T);
nt = getnt(T);
dt = getdt(T);

[b,A,B] = init_resolution(T,b,A,B);

switch lower(eulertype)
    case 'explicit'
        if isuniform(T) && ~istime(A)
            M = A/dt(1);
        else
            M = [];
        end
    case 'implicit'
        if isuniform(T) && ~istime(A) && ~istime(B)
            M = A/dt(1)+B;
        else
            M = [];
        end
end

if nargin>=6 && ~isempty(solver)
    solver = fcnchk(solver);
else
    if ~isempty(M) && isa(M,'double')
        [LM,UM] = lu(M);
        solver = @(A,b) solve(UM,(LM\b));
    else
        solver = @solve;
    end
end

n = size(A,1);
ut = cell(1,length(T));
vt = cell(1,length(T));
if nargin<5 || isempty(u0)
    ut{1} = zeros(n,1);
else
    ut{1} = u0;
end

if display_
    fprintf('Euler %s : ',eulertype);
end
switch lower(eulertype)
    case 'explicit'
        for i=1:nt
            if display_
                pourcentage(i,nt)
            end
            
            if istime(A)
                Ai = getmatrixatstep(A,i);
            else
                Ai = A;
            end
            if istime(B)
                Bi = getmatrixatstep(B,i);
            else
                Bi = B;
            end
            if ~isuniform(T) || istime(A)
                M = Ai/dt(i);
            end
            bi = b{i} + Ai*ut{i}/dt(i) - Bi*ut{i};
            [ut{i+1},result.solver{i}] = solver(M,bi);
            vt{i} = (ut{i+1} - ut{i})/dt(i);
            % vt{i} = solver(Ai,b{i}-Bi*ut{i});
        end
        vt{nt+1} = zeros(n,1);
        
    case 'implicit'
        vt{1} = zeros(n,1);
        for i=1:nt
            if display_
                pourcentage(i,nt)
            end
            
            if istime(A)
                Ai = getmatrixatstep(A,i);
            else
                Ai = A;
            end
            if istime(B)
                Bi = getmatrixatstep(B,i);
            else
                Bi = B;
            end
            if ~isuniform(T) || istime(A) || istime(B)
                M = Ai/dt(i)+Bi;
            end
            bi = b{i+1} + Ai*ut{i}/dt(i);
            [ut{i+1},result.solver{i}] = solver(M,bi);
            vt{i+1} = (ut{i+1} - ut{i})/dt(i);
            % vt{i+1} = solver(Ai,b{i+1}-Bi*ut{i+1});
        end
end

result.totaltime = toc(tStart);
if display_
    fprintf('Elapsed time = %f\n',result.totaltime)
end

try
    ut = horzcat(ut{:});
    vt = horzcat(vt{:});
catch
    warning('cell array')
end
ut = TIMEMATRIX(ut,T,[n,1]);
vt = TIMEMATRIX(vt,T,[n,1]);

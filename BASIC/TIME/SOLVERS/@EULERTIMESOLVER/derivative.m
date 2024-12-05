function vt = derivative(L,ut,b,A,B,solver)
% function vt = derivative(L,ut,b,A,B,solver)
% Computes the first order derivative v of u 
% by solving the differential equation Au'+Bu=b
% L: EULERTIMESOLVER
% b: right-hand side
% A,B: double or random matrices
% solver: function for solving linear systems
%         solver(A,b) solves Au=b
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'

eulertype = getparam(L,'eulertype');
T = gettimemodel(L);
nt = getnt(T);
dt = getdt(T);

% if nargin>=6 && ~isempty(solver)
%     solver = fcnchk(solver);
% else
%     if (~istime(A) && isa(A,'double')) || strcmpi(eulertype,'explicit')
%         [LA,UA] = lu(A);
%         solver = @(A,b) solve(UA,(LA\b));
%     else
%         solver = @solve;
%     end
% end

% n = size(A,1);
n = size(ut,1);
vt = cell(1,length(T));
if strcmpi(eulertype,'implicit')
    vt{1} = zeros(n,1);
end

for i=1:nt
%     if istime(A)
%         Ai = getmatrixatstep(A,i);
%     else
%         Ai = A;
%     end
%     if istime(B)
%         Bi = getmatrixatstep(B,i);
%     else
%         Bi = B;
%     end
    
    switch lower(eulertype)
        case 'explicit'
            vt{i} = (ut{i+1} - ut{i})/dt(i);
            % vt{i} = solver(Ai,b{i}-Bi*ut{i});
        case 'implicit'
            vt{i+1} = (ut{i+1} - ut{i})/dt(i);
            % vt{i+1} = solver(Ai,b{i+1}-Bi*ut{i+1});
    end
end

if strcmpi(eulertype,'explicit')
    vt{nt+1} = zeros(n,1);
end

try
    vt = horzcat(vt{:});
catch
    warning('cell array')
end
vt = TIMEMATRIX(vt,T,[n,1]);

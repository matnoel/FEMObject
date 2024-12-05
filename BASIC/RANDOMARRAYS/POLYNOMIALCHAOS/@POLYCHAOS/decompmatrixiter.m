function varargout = decompmatrixiter(apc,ng,h,fun,varargin)
% function u = decompmatrixiter(PC,ng,[],fun,varargin)
% Decomposition of a function of random variables on polynomial chaos basis
% by L2-projection
% PC : POLYCHOAS or PCMODEL or FESTOMODEL or PCMATRIX or PCARRAY
% ng : number of Gauss points, getorder(PC)+2 by default
% fun : function_handle of random variables RANDVARS(PC) associated to PC
% varargin : arguments of function fun
% call yk = fun(x(k,:),varargin{:}) for each sample vector x(k,:)
% x(:,i) is the vector containing all the evaluations of random variable i
% x(k,:) is the vector containing the k-th evaluation of all random variables
% yk is the response vector contaning the values returned by function fun applied to sample vector x(k,:)
% if PC is a POLYCHAOS, x is the sample matrix containing the evaluations of random variables associated to POLYCHAOS
%          a PCMODEL or a FESTOMODEL, x is the sample matrix containing the evaluations of random variables associated to PCMODEL
%          a PCMATRIX or a PCARRAY, x is the sample matrix containing the evaluations of random variables associated to PCMATRIX
%
% function u = decompmatrixiter(PC,ng,h,fun,varargin)
% h : RANDPOLYS used to compute the Gauss points

if ischarin('cell',varargin)
    varargin = delonlycharin('cell',varargin);
    varargout = decompmatrixiter_cell(apc,ng,h,fun,varargin{:});
    return
end

if ischarin('display',varargin)
    display_ = 1;
    varargin = delonlycharin('display',varargin);
else
    display_ = 0;
end

n = nargout;
fun = fcnchk(fun);
PC = getPC(apc);

if isempty(ng)
    ng = getorder(PC)+2;
elseif isa(ng,'function_handle')
    ng = ng(PC);
end

if isempty(h)
    gauss = calc_gausspoints(PC,ng);
    xpc = gauss.coord;
else
    h = RANDPOLYS(h);
    rv = RANDVARS(h);
    RV = RANDVARS(PC);
    gauss = calc_gausspoints(h,ng);
    xpc = transfer(rv,RV,gauss.coord);
end

N = gauss.nbgauss;
if display_
    fprintf('\nNumber of Gauss points = %d\n',N);
end

% x = transfer(apc,xpc);
% t0 = clock;
% for k=1:N
%     if display_
%         pourcentage(k,N)
%     end
%     Ak = polyval(PC,xpc(k,:));
%     wk = gauss.w(k);
%     yk = cell(1,n);
%     [yk{:}] = fun(x(k,:),varargin{:});
%     for i=1:n
%         if isa(yk{i},'cell')
%             for l=1:length(yk{i})
%                 uk = wk*yk{i}{l}(:)*Ak;
%                 if k==1
%                     u{i}{l} = uk;
%                     s{i}{l} = size(yk{i}{l});
%                 else
%                     u{i}{l} = u{i}{l} + uk;
%                 end
%             end
%         else
%             uk = wk*yk{i}(:)*Ak;
%             if k==1
%                 u{i} = uk;
%                 s{i} = size(yk{i});
%             else
%                 u{i} = u{i} + uk;
%             end
%         end
%         
%     end
%     if display_ && k==3
%         fprintf('\nNumber of Gauss points = %d -- Estimated time : %4.1f s\n',N,etime(clock,t0)*N/k);
%     end
% end

A = polyval(PC,xpc);
W = spdiags(gauss.w(:),0,N,N);
x = transfer(apc,xpc);
[y,s] = calc_response_iter(fun,x,n,varargin{:});

u = cell(1,n);
for i=1:n
    if isa(y{i},'cell')
        for l=1:length(y{i})
            u{i}{l} = y{i}{l}*W*A;
        end
    else
        u{i} = y{i}*W*A;
    end
end

for i=1:n
    if isa(u{i},'cell')
        for l=1:length(u{i})
            u{i}{l} = PCMATRIX(u{i}{l},s{i}{l},PC);
        end
    else
        u{i} = PCMATRIX(u{i},s{i},PC);
    end
end

varargout = u;

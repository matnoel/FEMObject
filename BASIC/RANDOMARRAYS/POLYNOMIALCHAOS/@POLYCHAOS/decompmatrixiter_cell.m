function u = decompmatrixiter_cell(apc,ng,h,fun,varargin)
% function u = decompmatrixiter_cell(PC,ng,[],fun,varargin)
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
% function u = decompmatrixiter_cell(PC,ng,h,fun,varargin)
% h : RANDPOLYS used to compute the Gauss points

if ischarin('display',varargin)
    display_ = 1;
    varargin = delonlycharin('display',varargin);
else
    display_ = 0;
end

fun = fcnchk(fun);
PC = getPC(apc);

if isempty(ng)
    ng = getorder(PC)+2;
end

if isempty(h)
    gauss = calc_gausspoints(PC,ng);
    x = gauss.coord;
else
    h = RANDPOLYS(h);
    rv = RANDVARS(h);
    RV = RANDVARS(PC);
    gauss = calc_gausspoints(h,ng);
    x = transfer(rv,RV,gauss.coord);
end

N = gauss.nbgauss;
if display_
    fprintf('\nNumber of Gauss points = %d\n',N);
end

xpc = x;
x = transfer(apc,xpc);

t0 = clock;
u = cell(1,length(PC));
for k=1:N
    if display_
        pourcentage(k,N)
    end
    Dk = polyval(PC,xpc(k,:));
    wk = gauss.w(k);
    yk = fun(x(k,:),varargin{:});
    s = size(yk);
    for i=1:length(PC)
        uk = wk*yk*Dk(i);
        if k==1
            u{i} = uk;
        else
            u{i} = u{i} + uk;
        end
    end
    if display_ && k==3
        fprintf('\nNumber of Gauss points = %d -- Estimated time : %4.1f s\n',N,etime(clock,t0)*N/k);
    end
end

u = PCMATRIX(u,s,PC);

% u = cell(1,length(PC));
% for i=1:length(PC)
%     u{i} = sparse((gauss.w.*D(:,i)'))*y;
%     u{i} = reshape(u{i},s);
% end
% u = PCCELL(u,POLYCHAOS(PC));


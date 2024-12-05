function varargout = decompmatrix(apc,ng,h,fun,varargin)
% function u = decompmatrix(PC,ng,[],fun,varargin)
% Decomposition of a function of random variables on polynomial chaos basis
% by L2-projection
% PC : POLYCHOAS or PCMODEL or FESTOMODEL or PCMATRIX or PCARRAY
% ng : number of Gauss points, getorder(PC)+2 by default
% fun : function_handle of random variables RANDVARS(PC) associated to PC
% varargin : arguments of function fun
% call y = fun(x,varargin{:}) for sample matrix x
% x(:,i) is the sample vector containing all the evaluations of random variable i
% x(k,:) is the sample vector containing the evaluation k of all random variables
% y is the response matrix contaning the values returned by function fun applied to sample matrix x
% if PC is a POLYCHAOS, x is the sample matrix containing the evaluations of random variables associated to POLYCHAOS
%          a PCMODEL or a FESTOMODEL, x is the sample matrix containing the evaluations of random variables associated to PCMODEL
%          a PCMATRIX or a PCARRAY, x is the sample matrix containing the evaluations of random variables associated to PCMATRIX
%
% function u = decompmatrix(PC,ng,h,fun,varargin)
% h : RANDPOLYS used to compute the Gauss points

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

A = polyval(PC,xpc);
W = spdiags(gauss.w(:),0,N,N);
x = transfer(apc,xpc);
[y,s] = calc_response(fun,x,n,varargin{:});

u = cell(1,n);
for i=1:n
    if isa(y{i},'cell')
        for l=1:length(y{i})
            u{i}{l} = y{i}{l}*W*A;
            u{i}{l} = PCMATRIX(u{i}{l},s{i}{l},PC);
        end
    else
        u{i} = y{i}*W*A;
        u{i} = PCMATRIX(u{i},s{i},PC);
    end
end

varargout = u;

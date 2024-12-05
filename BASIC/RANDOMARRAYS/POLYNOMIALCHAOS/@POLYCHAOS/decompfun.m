function u = decompfun(apc,ng,h,fun,varargin)
% function u = decompfun(PC,ng,[],fun,varargin)
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
% function u = decompfun(PC,ng,h,fun,varargin)
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
elseif isa(ng,'function_handle')
    ng = ng(PC);
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

D = polyval(PC,x);
x = transfer(apc,x);
y = fun(x,varargin{:});

if numel(y)~=N
    error('pour decomposer un vecteur ou une matrice , utiliser decompmatrix')
end

y = double(y);
W = spdiags(gauss.w(:),0,N,N);
u = y(:)'*W*D;
u = PCMATRIX(u,[1,1],PC);


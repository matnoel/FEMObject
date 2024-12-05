function varargout = decompmatrix_collocation(apc,fun,varargin)
% function u = decompmatrix_collocation(PC,fun,varargin)
% Decomposition of a function of random variables on polynomial chaos basis
% by Lagrange interpolation (or collocation)
% PC : POLYCHOAS or PCMODEL or FESTOMODEL or PCMATRIX or PCARRAY
% fun : function_handle of random variables RANDVARS(PC) associated to PC
% varargin : arguments of function fun
% call y = fun(x,varargin{:}) for sample matrix x
% x(:,i) is the sample vector containing all the evaluations of random variable i
% x(k,:) is the sample vector containing the evaluation k of all random variables
% y is the response matrix contaning the values returned by function fun applied to sample matrix x
% if PC is a POLYCHAOS, x is the sample matrix containing the evaluations of random variables associated to POLYCHAOS
%          a FESTOMODEL or a PCMODEL, x is the sample matrix containing the evaluations of random variables associated to PCMODEL
%          a PCMATRIX or a PCARRAY, x is the sample matrix containing the evaluations of random variables associated to PCMATRIX

if ischarin('display',varargin)
    display_ = 1;
    varargin = delonlycharin('display',varargin);
else
    display_ = 0;
end

n = nargout;
fun = fcnchk(fun);
PC = getPC(apc);

xpc = calc_collocationpoints(PC);

N = size(x,1);
if display_
    fprintf('\nNumber of collocation points = %d\n',N);
end

x = transfer(apc,xpc);
[u,s] = calc_response(fun,x,n,varargin{:});

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

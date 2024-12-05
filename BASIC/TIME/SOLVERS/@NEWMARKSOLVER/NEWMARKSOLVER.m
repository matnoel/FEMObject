function L = NEWMARKSOLVER(T,varargin)
% function L = NEWMARKSOLVER(T,'alpha',alpha)
% T : TIMEMODEL
% alpha : parameter (0.05 by default)
% gamma : parameter (1/2+alpha by default)
% beta : parameter ((1+alpha)^2/4 by default)
% alpha=0, gamma=0,   beta=0    : explicit method (unstable)
% alpha=0, gamma=1/2, beta=0    : explicit central-difference method (conditionally stable)
% alpha=0, gamma=1/2, beta=1/12 : implicit Fox and Godwin method (conditionally stable)
% alpha=0, gamma=1/2, beta=1/6  : implicit linear acceleration method (conditionnally stable)
% alpha=0, gamma=1/2, beta=1/4  : implicit mean acceleration method (unconditionnally stable)
% alpha>0, gamma=1/2+alpha, beta=(1+alpha)^2/4 : implicit modified acceleration method (unconditionnally stable)
% If beta=0, the Newmark method is explicit. 
% If gamma<1/2, the Newmark method is unstable.
% If gamma>=1/2 and 2*beta>=gamma, the Newmark method is unconditionnaly stable (no restriction on the time step dt).
% If gamma>=1/2 and 2*beta<gamma, the Newmark method is conditionnaly stable.

if nargin==0
    L = struct();
    L = class(L,'NEWMARKSOLVER',PARAMETERS(),TIMEMODEL());
else
    param = PARAMETERS(varargin{:});
    param = setdefaultparam(param,'alpha',0.05);
    alpha = getparam(param,'alpha');
    param = setdefaultparam(param,'gamma',1/2+alpha);
    param = setdefaultparam(param,'beta',(1+alpha)^2/4);
    param = setdefaultparam(param,'display',false);
    
    T = setapproxparam(T,'p',1);
    T = setapproxparam(T,'type','default');
    
    L = struct();
    L = class(L,'NEWMARKSOLVER',param,T);
end

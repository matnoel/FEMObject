function L = DGTIMESOLVER(T,p,varargin)
% function L = DGTIMESOLVER(T,p,varargin)
% T : TIMEMODEL
% p : polynomial degree

if nargin==0
    L = struct();
    L = class(L,'DGTIMESOLVER',PARAMETERS(),TIMEMODEL());
else
    if nargin==1 || ~isa(p,'double')
        error('rentrer le degre de l''approximation polynomiale')
    end
    param = PARAMETERS(varargin{:});
    param = setdefaultparam(param,'outputsplit',true);
    param = setdefaultparam(param,'display',false);
    param = setdefaultparam(param,'lu',true);
    
    checkparamvalue(param,'outputsplit',true,false);
    checkparamvalue(param,'display',true,false);
    checkparamvalue(param,'lu',true,false);
    
    T = setapproxparam(T,'p',p);
    T = setapproxparam(T,'type','DG');
    
    L = struct();
    L = class(L,'DGTIMESOLVER',param,T);
end

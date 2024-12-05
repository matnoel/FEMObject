function L = EULERTIMESOLVER(T,varargin)
% function L = EULERTIMESOLVER(T,'eulertype',eulertype)
% T : TIMEMODEL
% eulertype : 'implicit' or 'explicit' (implicit by default)
%
% EULERTIMESOLVER(T,'display',true)
% display the progression of the time resolution (false by default)

if nargin==0
    L = struct();
    L = class(L,'EULERTIMESOLVER',PARAMETERS(),TIMEMODEL());
else
    param = PARAMETERS(varargin{:});
    param = setdefaultparam(param,'eulertype','implicit');
    param = setdefaultparam(param,'display',false);
    
    checkparamvalue(param,'eulertype','implicit','explicit');
    checkparamvalue(param,'display',true,false);
    
    T = setapproxparam(T,'p',1);
    T = setapproxparam(T,'type','default');
    
    L = struct();
    L = class(L,'EULERTIMESOLVER',param,T);
end

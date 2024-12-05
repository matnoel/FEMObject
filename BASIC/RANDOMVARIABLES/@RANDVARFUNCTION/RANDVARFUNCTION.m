function rv = RANDVARFUNCTION(fun,varargin)
% function rv = RANDVARFUNCTION(fun,varargin)
% varargin : arguments de la function fun
% fun inline ou function_handle avec autant d'arguments que de RANDVARS

rv.param = cell(0,1);
rv.randomparam = [];
for i=1:length(varargin)
    rv.param{i}=varargin{i};   
    if isa(varargin{i},'RANDVAR')
        rv.randomparam=[rv.randomparam,i];
    elseif isa(varargin{i},'RANDVARS')
        error('les arguments de la function ne doivent pas etre des RANDVARS');    
    end
end
rv.RV = RANDVARS(rv.param{rv.randomparam});
if length(rv.RV)==0
    warning('aucun argument aleatoire');
    rv = fun(rv.param{:});
    return
end
if ~(isa(fun,'function_handle') | isa(fun,'inline'))
    error('premier argument doit etre une function_handle ou inline')    
end

rv.fun = fun;

rv = class(rv,'RANDVARFUNCTION');


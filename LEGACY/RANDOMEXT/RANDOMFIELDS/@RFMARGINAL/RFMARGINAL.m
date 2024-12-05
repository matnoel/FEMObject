function ma = RFMARGINAL(rv,varargin)
% function ma = RFMARGINAL(rv)
% rv : random variable
% 
% function ma = RFMARGINAL(rv,paramname1,param1fun,paramname2,param2fun,...)
% to enter a function for the parameter of the random variable
% paramnamei : name of the random variable parameter
% param1fun : inline function ou char expression of the function
%   functions associated with different arguments must have the same
%   parameters and order of paramaters (for the eval function)
% be careful of the char expression : arguments are seached with symvar
% so if more than one argument, use an inline function for input

if nargin==0
    ma=struct('RV',RANDVAR(),'param',struct());
    ma = class(ma,'RFMARGINAL');       
else
    ma.RV = rv;
    param = getparam(rv);
    ma.param=param;
    if nargin>1 
        paramnames = fieldnames(param);

        for k=1:length(paramnames)
            val = getcharin(paramnames{k},varargin);
            if ~isempty(val)
                [rep,pos]=ischarin(paramnames{k},varargin);
                if isa(val,'char')
                    var = symvar(val);
                    val = inline(val,var{:});
                end
                ma.param=setfield(ma.param,paramnames{k},val);
            end
        end
    end
    ma = class(ma,'RFMARGINAL');       
end


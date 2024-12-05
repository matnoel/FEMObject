function plot_lssol(elem,node,q,ls,varargin)
% function plot_lssol(elem,node,q,ls,varargin)

if isempty(getlsnumber(elem))
    nature = 'domain';
    ls = [];
else
    ls = getlevelset(ls,getlsnumber(elem));
    nature = getnature(ls);
end

switch nature
    case 'crack'
        if ~isenrich(elem)
            plot_sol(elem,node,q,varargin{:});
        else
            plot_lscracksol(elem,node,q,ls,varargin{:});
        end
    case 'material'
        plot_lsmaterialsol(elem,node,q,ls,varargin{:});
    case 'domain'
        plot_lsdomainsol(elem,node,q,ls,varargin{:});
end

end


function plot_enrichfun(elem,node,ls,lsnum,varargin)

ls = getlevelset(ls,lsnum);
nature = getnature(ls);

switch nature
    case 'material'
        if isenrich(ls)
 plot_enrichfunmaterial(elem,node,ls,varargin{:}); 
        end
    otherwise
 
end

end


function ke = rigilspc(elem,node,PC,ls,varargin)


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
          ke = rigipc(elem,node,PC,varargin{:});
        else
          ke = rigilspccrack(elem,node,PC,ls,varargin{:});
        end
    case 'material'
        ke = rigilspcmaterial(elem,node,PC,ls,varargin{:});
    case 'domain'       
        ke = rigilspcdomain(elem,node,PC,ls,varargin{:});       
    otherwise 
        error('pas defini')
end

function ke = rigils(elem,node,ls,varargin)

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
          ke = rigi(elem,node,varargin{:});
        else
          ke = rigilscrack(elem,node,ls,varargin{:});
        end
    case 'material'
        ke = rigilsmaterial(elem,node,ls,varargin{:});
    case 'domain'
        ke = rigilsdomain(elem,node,ls,varargin{:});
    otherwise 
        error('pas defini')
end

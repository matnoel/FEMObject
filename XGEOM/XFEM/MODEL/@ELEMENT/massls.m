function me = massls(elem,node,ls,varargin)

if isempty(getlsnumber(elem))
    nature = 'domain';
    ls = [];  
else
    ls = getlevelset(ls,getlsnumber(elem));
    nature = getnature(ls);
end

switch nature

    case 'domain'
        me = masslsdomain(elem,node,ls,varargin{:});
    otherwise 
        error('pas defini')
end

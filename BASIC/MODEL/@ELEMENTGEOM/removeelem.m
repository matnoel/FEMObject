function elem=removeelem(elem,num,varargin)
% function elem=removeelem(elem,num)
% num est la liste des numeros des elements (numerotation locale au groupe d'elements)
% function elem=removeelem(elem,num,'global')
% num est la liste des numeros des elements (numerotation globale)

if nargin>2 &  strcmp(varargin{1},'global')
    ia = find(ismember(getnumber(elem),num));
    varargin{1}='local';
else
    ia = num;
end

elem = getelem(elem,setdiff(1:getnbelem(elem),ia),varargin{:});


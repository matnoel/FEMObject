function elem = STOKES(node,numelem,connec,varargin)
% function elem = STOKES(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements

if nargin==0
    elemp = TRI3();
else
    elemp = TRI3(node,numelem,connec,varargin{:});
end
elem = struct();
elem = class(elem,'STOKES',elemp);
elem = setlocal(elem);

function elem = COQ4(node,numelem,connec,varargin)
% function elem = COQ4(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements
% function elem = COQ4(node,numelem,connec,...'fullintegration') 

if nargin==0
    elemp = QUA4();
else
    elemp = QUA4(node,numelem,connec,varargin{:});
end
elem = struct();
elem = class(elem,'COQ4',elemp);
elem = setlocal(elem);
elem = setparam(elem,'fullintegration',ischarin('fullintegration',varargin));

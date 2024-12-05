function elem = DKQ(node,numelem,connec,varargin)
% function elem = DKQ(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements

if nargin==0
    elemp = QUA4();
else
    elemp = QUA4(node,numelem,connec,varargin{:});
end
elem = struct();
elem = class(elem,'DKQ',elemp);
elem = setlocal(elem);

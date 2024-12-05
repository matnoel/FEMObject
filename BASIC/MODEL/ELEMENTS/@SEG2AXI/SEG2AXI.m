function elem = SEG2AXI(node,numelem,connec,varargin)
% function elem = SEG2AXI(node,numelem,connec,'material',mat,'option',option)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements

if nargin==0
    elemp = SEG2();
else
    elemp = SEG2(node,numelem,connec,varargin{:});
end
elem = struct();
elem = class(elem,'SEG2AXI',elemp);

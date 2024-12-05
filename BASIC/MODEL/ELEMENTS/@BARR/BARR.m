function elem = BARR(node,numelem,connec,varargin)
% function elem = BARR(node,numelem,connec,'material',mat)
% node : objet de type NODE contenant les noeuds de l'element
% mat : materiau
% connec : table de connectivite
% numelem : numero des elements

if nargin==0
    elemp = SEG2();
else
    elemp = SEG2(node,numelem,connec,varargin{:});
end
elem = struct();
elem = class(elem,'BARR',elemp);

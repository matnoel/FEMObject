function elem = SEG2DK(varargin)
% function elem = SEG2DK(node,numelem,connec,'material',mat,'option',option,syscoord,'parent',parent)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements

if nargin==0
    elemp = SEG2();
elseif nargin==1 && isa(varargin{1},'SEG2')
    elemp = varargin{1};
else
    elemp = SEG2(varargin{:});
end
elem = struct();
elem = class(elem,'SEG2DK',elemp);
elem = setlocal(elem);

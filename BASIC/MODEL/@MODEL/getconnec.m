function connec = getconnec(M,type)
% function connec = getconnec(M,type)
% type = 'node2node' ou 'node2elem' ou 'elem2node' ou 'elem2elem'
% renvoi la matrice correspondante dans la structure M.connec
%
% function connec = getconnec(M)
% renvoi la structure connec
%
% See also MODEL/calc_connec


connec = M.connec;
if nargin==2
    connec = eval(['connec.' type]);
end

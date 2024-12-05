function ddl = createddlvect(S,name,type)
% function ddl = createddlvect(S,name,type)
% name : nom du vecteur
% S : systeme de coordonnes
% type : 'ROTA' ou 'TRANS'

if nargin<=2
    type = 'TRANS';
end
switch type
    case 'TRANS'
        ddl = {[name 'X'],[name 'Y']};
    case 'ROTA'
        ddl = {[name 'Z']};
end

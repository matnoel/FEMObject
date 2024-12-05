function A=noenrichmatrix(BC,A,choix)

if nargin==2
    A = A(BC.ddlnoenrich,BC.ddlnoenrich);
elseif strcmp(choix,'free')
    rep = getddlnoenrich(BC,'free')
    A = A(rep,rep);
elseif strcmp(choix,'bloque')
    rep = getddlnoenrich(BC,'bloque')
    A = A(rep,rep);
end

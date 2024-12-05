function d=getddlnoenrich(B,choix)

if nargin==1
d=B.ddlnoenrich;
elseif strcmp(choix,'free')
[a,d]=ismember(B.ddlnoenrich,B.ddlfree);
d=nonzeros(d);
elseif strcmp(choix,'bloque')
[a,d]=ismember(B.ddlnoenrich,B.ddlbloque);
d=nonzeros(d);
end

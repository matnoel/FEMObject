function d=getddlenrich(B,choix)

if nargin==1
d=B.ddlenrich;
elseif strcmp(choix,'free')
[a,d]=ismember(B.ddlenrich,B.ddlfree);
d=nonzeros(d);
elseif strcmp(choix,'bloque')
[a,d]=ismember(B.ddlenrich,B.ddlbloque);
d=nonzeros(d);
end
    
    
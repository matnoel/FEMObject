function A=enrichmatrix(BC,A,choix)

if nargin==2
A = A(BC.ddlenrich,BC.ddlenrich);
elseif strcmp(choix,'free')
rep = getddlenrich(BC,'free')
A = A(rep,rep);
elseif strcmp(choix,'bloque')
rep = getddlenrich(BC,'bloque')
A = A(rep,rep);
end

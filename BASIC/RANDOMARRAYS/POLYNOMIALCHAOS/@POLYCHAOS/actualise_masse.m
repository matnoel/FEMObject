function PC=actualise_masse(PC,PC2)


if nargin==1
if ~isa(PC.masse,'MULTIMATRIX') || length(PC.masse)~=length(PC)
try
PC.masse = PC.masse{1:length(PC)}(1:length(PC),1:length(PC));
catch
%warning('calcul de la matrice masse stochastique')
PC = calc_masse(PC);
end
end
else
if ~isa(PC.masse,'MULTIMATRIX') || length(PC.masse)~=length(PC) || size(PC.masse,1)~=length(PC2)    
try    
if polycmp(PC,PC2)
PC.masse = PC.masse{1:length(PC)}(1:length(PC2),1:length(PC2));
else
[rep,ia] = isin(PC,PC2);
if rep
PC.masse = PC2.masse{ia}(1:length(PC2),1:length(PC2)) ;
else
[rep,ia] = isin(PC2,PC);
PC.masse = PC.masse(ia,ia) ;
end
end
catch
%warning('calcul de la matrice masse stochastique')
PC = calc_masse(PC,PC2);    
end
end
end
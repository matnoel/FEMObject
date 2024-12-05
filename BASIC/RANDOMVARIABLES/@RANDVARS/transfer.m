function bx = transfer(a,b,ax)
% function bx = transfer(a,b,ax)
% a et b : RANDVARS
% les dimensions stochastiques de b doivent etre contenues dans celles de a
%     -> utilisation de ismember(b,a)
% variables supposees independantes
% nombre de variables a = Ma
% nombre de variables b = Mb (b peut etre une RANDVAR)
%
% ax : n-by-Ma double ou 1-by-Ma cell of n-by-1 double
% bx : n-by-Mb double
%
% Pour deux variables ai et bj de meme dimension stochastique
% bxj = F_bj^-1(F_ai(axi))
%
% See also RANDVAR/transfer, RANDVAR/randomeval, RANDVARS/randomeval

[ok,rep] = ismember(b,a);
if ~all(ok)
    error(['les dimensions stochastiques de b doivent etre incluses dans celles de a'])
end

if isa(ax,'cell')
    ax = [ax{:}];
end

if size(ax,2)~=length(a)
    error(['ax doit avoir le meme nombre de colonnes (ou de cellules) que le nombre de variables aleatoires a'])
end

a = a.RV ;
if isa(b,'RANDVAR')
    b = {b};
else
    b = b.RV;
end

bx = zeros(size(ax,1),length(b));

for i=1:length(rep)
    bx(:,i) = transfer(a{rep(i)},b{i},ax(:,rep(i)));
end

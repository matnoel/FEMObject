function hval = dpolyval(h,liste,x)

% hval=dpolyval(h,liste,x)
%
% Evaluate a the derivative of a POLYFELAGRANGE h_n(x)
%
% x : vecteur (ou matrice) des points ou les polynomes sont evaluees
% liste : liste des indices des polynomes a evaluer
% hval = matrice des valeurs . hval(i,j) = valeur de h_liste(j) en x(i)
%
% les poynomes de lagrange sont tries par ordre croissant de la position
% des points d'interpolation


switch min(size(x))
    case 0
hval=sparse(0,length(liste));     
    case 1
param = get(h,'param');
n=param.n;
p=param.p;
m=param.m;
subpoints = param.subpoints;
reppoints = param.reppoints;

hval=zeros(numel(x),length(liste));
x=x(:);

for i=1:n
    x1 = param.I(i,1);
    x2 = param.I(i,2);
    dx=x2-x1;
 
    if i<n
    rep = find(x>=x1 & x<x2);
    else
    rep = find(x>=x1 & x<=x2);
    end
    [rep1,rep2] = ismember(liste+1,reppoints{i});
    rep1 = find(rep1);
    rep2 = rep2(rep1);
 hval(rep,rep1)=dpolyval(POLYLAGRANGE(subpoints{i}),rep2-1,x(rep));
 
 end 
hval = sparse(hval) ;
    otherwise
hval = polyval(h,liste,x(:));
hval = reshape(full(hval),[size(x,1),size(x,2),length(liste)]);

end
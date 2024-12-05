function c = polycoeff(h,liste)

% c=polycoeff(h,liste)
%
% Coefficient of Laguerre polynomials L^a_n
% polynomes orthonormés sur R+ pour la mesure x^(a-1)*exp(-x)/Gamma(a)
%
% entree : liste = liste des indices des polynomes
% sortie : c = matrice des coefficients c(i+1,j+1) = coefficient en x^j de L^a_(liste(i))

param=get(h,'param');
a=param.a;

i=max(liste);
c=zeros(i+1,i+1);

c(1,1)=1;
if i>=1
c(2,2)=-1;c(2,1)=a;
if i>=2
for n=2:i
c(n+1,2:n+1)=-1/n*c(n,1:n);
c(n+1,1:n)=c(n+1,1:n)+(2+(a-2)/n)*c(n,1:n);
c(n+1,1:n-1)=c(n+1,1:n-1)-(1+(a-2)/n)*c(n-1,1:n-1);
end
end

end

listeu = unique(liste);
for k=listeu(:)'
c(k+1,:)=c(k+1,:)/sqrt(gamma(k+a)/gamma(a)/factorial(k)); % poly normalise
end
c=c(liste+1,:);

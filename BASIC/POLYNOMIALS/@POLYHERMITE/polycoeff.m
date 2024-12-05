function c = polycoeff(h,liste)

% c=polycoeff(h,liste)
%
% Coefficient of Hermite polynomials He_n
% polynomes orthonorm�s sur R pour la mesure 1/sqrt(2*pi)*exp(-x^2/2)
%
% liste : liste des indices des polynomes
% c : matrice des coefficients c(i+1,j+1) = coefficient en x^j de He_(liste(i))


i=max(liste);
c=zeros(i+1,i+1);

c(1,1)=1;
if i>=1
c(2,2)=1;
if i>=2
for n=2:i
c(n+1,2:n+1)=c(n,1:n);
c(n+1,1:n-1)=c(n+1,1:n-1)-(n-1)*c(n-1,1:n-1);
end
end

end

listeu = unique(liste);
for k=listeu(:)'
c(k+1,:)=c(k+1,:)/sqrt(factorial(k)); % poly normalise
end
c=c(liste+1,:);


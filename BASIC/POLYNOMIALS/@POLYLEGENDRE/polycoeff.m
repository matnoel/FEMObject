function c = polycoeff(h,liste)

% c=legendrecoeff(h,liste)
%
% Coefficient of Legendre polynomials Le_n
% polynomes orthonormés sur [-1,1] pour la mesure 1/2
%
% entree : liste : liste des indices des polynomes
% sortie : c : matrice des coefficients c(i+1,j+1) = coefficient en x^j de Le_(liste(i))


i=max(liste);
c=zeros(i+1,i+1);

c(1,1)=1;
if i>=1
c(2,2)=1;
%c(2,2)=sqrt(3);
if i>=2
for n=2:i
c(n+1,2:n+1)=(2-1/n)*c(n,1:n);
c(n+1,1:n-1)=c(n+1,1:n-1)-(1-1/n)*c(n-1,1:n-1);
%c(n+1,2:n+1)=sqrt(4-1/n^2)*c(n,1:n);
%c(n+1,1:n-1)=c(n+1,1:n-1)-sqrt((2*n+1)/(2*n-3))*(1-1/n)*c(n-1,1:n-1);
end
end

end

listeu = unique(liste);
for k=listeu(:)'
c(k+1,:)=c(k+1,:)*sqrt(2*k+1); % poly normalise
end
c=c(liste+1,:);

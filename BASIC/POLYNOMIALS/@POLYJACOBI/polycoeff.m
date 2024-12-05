function c = polycoeff(h,liste)

% c=polycoeff(h,liste)
%
% Coefficient of jacobi polynomials P(a,b)_n(x)
% polynomes orthonormés sur [-1,1] pour la mesure (1-x)^(b-1)*(1+x)^(a-1)/Beta(a,b)/2^(a+b-1)
%
% entree : liste = liste des indices des polynomes
% sortie : c = matrice des coefficients c(i+1,j+1) = coefficient en x^j de P(a,b)_(liste(i))
param=get(h,'param');
a=param.a-1;
b=param.b-1;

i=max(liste);
c=zeros(i+1,i+1);

c(1,1)=1;
if i>=1
c(2,2)=(2+a+b)/2;c(2,1)=(b-a)/2;
if i>=2
for n=2:i
d1=(2*n+a+b-1)*(2*n+a+b)/(2*n)/(n+a+b);
d2=(a^2-b^2)*(2*n+a+b-1)/(2*n+a+b-2)/(2*n)/(n+a+b);
d3=(n+b-1)*(n+a-1)*(2*n+a+b)/(2*n+a+b-2)/(n)/(n+a+b);
% P(n)=(d1*x-d2)P(n-1) - d3P(n-2)
c(n+1,2:n+1)=d1*c(n,1:n);
c(n+1,1:n)=c(n+1,1:n)-d2*c(n,1:n);
c(n+1,1:n-1)=c(n+1,1:n-1)-d3*c(n-1,1:n-1);
end
end

end

listeu = unique(liste);
for k=listeu(:)'
betak=1/(2*k+a+b+1)/factorial(k)*gamma(k+a+1)*gamma(k+b+1)/gamma(k+a+b+1)/beta(a+1,b+1);
c(k+1,:)=c(k+1,:)/sqrt(betak); % poly normalise
end
c=c(liste+1,:);

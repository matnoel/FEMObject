function c = polycoeff(h,liste)
% c=polycoeff(h,liste)
%
% Coefficient of orthonormal polynomials P_n(x) with respect to inner product
% int(P_n P_m p),  given
% int(P_n x^k p) for all k, where p is a probability density function
% entree : liste = liste des indices des polynomes
% sortie : c = matrice des coefficients c(i+1,j+1) = coefficient en x^j de P(a,b)_(liste(i))

%param=get(h,'param');

i=max(liste);
c=zeros(i+1,i+1);

intxn = calc_intxn(h,0:2*i);
IJ = repmat((0:i)',1,i+1)+repmat((0:i),i+1,1);

c(1,1)=1;
for n=1:i
    c(n+1,n+1) = 1 ;
    for k=0:n-1
        c(n+1,1:k+1) = c(n+1,1:k+1) - (intxn((n+1):(n+k+1))*c(k+1,1:k+1)').*c(k+1,1:k+1);
    end
    a = c(n+1,1:n+1)*intxn(1+IJ(1:n+1,1:n+1))*c(n+1,1:n+1)';
    c(n+1,:)=c(n+1,:)/sqrt(abs(a));
end

c=c(liste+1,:);

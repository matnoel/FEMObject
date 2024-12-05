function x = eicdf(Xs,y)
% function x = eicdf(Xs,y)
% Xs : echantillons d'une variable aleatoire X
% calcul de x = F^-1(y), y dans [0,1]
% o� F est la fonction de distribution de X


[F,X] = ecdf(Xs);

if any(y>1 | y<0)
    error('la icdf n''est definie que sur [0,1]')
end
x=zeros(size(y));
for i=1:length(y)
    r1 = max(find(F<y(i) | F==0));
    r2 = min(find(F>=y(i)));
    F1 = F(r1); 
    F2 = F(r2);
    x1 = X(r1);
    x2 = X(r2);
    x(i)=x2;
    x(i) = x1 + (y(i)-F1)/(F2-F1)*(x2-x1);   
end

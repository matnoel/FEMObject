function dhval = dpolyval(h,liste,x)

% dhval=dpolyval(h,liste,x)
%
% Evaluate the derivative of RANDOM POLYNOMIAL dh_n(x)/dx
%
% x : vecteur (ou matrice) des points ou les polynomes sont evaluees
% liste : liste des indices des polynomes a evaluer
% dhval = matrice des valeurs . dhval(i,j+1) = valeur de dh_liste(j) en x(i)

L = fliplr(dpolycoeff(h,liste));

switch min(size(x))
    case 1
dhval=zeros(numel(x),length(liste));
for i=1:length(liste)
dhval(:,i) = polyval(L(i,:),x(:));
end
    otherwise
dhval=zeros(size(x,1),size(x,2),length(liste));
for i=1:length(liste)
dhval(:,:,i) = polyval(L(i,:),x);
end
        
end
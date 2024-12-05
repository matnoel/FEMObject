function hval = polyval(h,liste,x)
% hval = polyval(h,liste,x)
%
% Evaluate a RANDOM POLYNOMIAL h_n(x)
%
% x : vecteur (ou matrice) des points ou les polynomes sont evaluees
% liste : liste des indices des polynomes a evaluer
% hval = matrice des valeurs . hval(i,j+1) = valeur de h_liste(j) en x(i)

L = fliplr(polycoeff(h,0:max(liste)));

switch min(size(x))
    case 0
        hval=sparse(0,length(liste));
    case 1
        hval=sparse(numel(x),size(L,1));
        for i=1:size(L,1)
            hval(:,i) = polyval(L(i,:),x(:));
        end
        
        hval = hval(:,liste+1);
    otherwise
        error('programmer en MULTIMATRIX')
        hval=zeros(size(x,1),size(x,2),length(liste));
        for i=1:length(liste)
            hval(:,:,i) = polyval(L(i,:),x);
        end
        
end
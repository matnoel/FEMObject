function hval = polyval(h,liste,x)
% hval = polyval(h,liste,x)
%
% Evaluate a RANDOM POLYNOMIAL h_n(x)
%
% x : vecteur (ou matrice) des points ou les polynomes sont evaluees
% liste : liste des indices des polynomes a evaluer
% hval = matrice des valeurs . hval(i,j) = valeur de h_liste(j) en x(i)

param = get(h,'param');
n=param.n;
% p=param.p;

switch min(size(x))
    case 0
        hval=sparse(0,length(liste));
    case 1
        hval=zeros(numel(x),length(liste));
        x=x(:);
        x=round(x*1e15)*1e-15;
        I = mod(liste,n)+1;
        P = (liste+1-I)/n ;
        
        for i=1:n
            x1 = param.I(i,1);
            x2 = param.I(i,2);
            dx=x2-x1;
            
            if i<n
                rep = find(x>=x1 & x<x2);
            else
                rep = find(x>=x1 & x<=x2);
            end
            repliste = find(I==i);
            
            
            xlocal = transfer(RVUNIFORM(x1,x2),RVUNIFORM(-1,1),x(rep));
            
            hval(rep,repliste)=polyval(POLYLEGENDRE(),P(repliste),xlocal)/sqrt(dx);
        end
        hval = sparse(hval) ;
    otherwise
        hval = polyval(h,liste,x(:));
        hval = reshape(full(hval),[size(x,1),size(x,2),length(liste)]);
        
end
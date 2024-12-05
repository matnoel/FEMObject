function phi = eval_ortho_matrix_param_inverse(n,p,A,B)
%  phi = eval_ortho_matrix_param(n,p,A)
%  evaluation des parametres d'une matrice orthogonale 
% A : matrice orthogonale n-by-p
%  phi : parametres 
% 
%  - cas n=1 : vecteur orthonorme dans R^p
%             les parametres sont les parametres angulaires d'une
%             hypersphere
%             phi = (t1,t2,...,tp-1) in [0,2pi]x[0,pi]x...x[0,pi]
%
%  - cas general variete de Stiefel S(p,n)
% 
%  - cas n=p : groupe orthogonal O+(n)
%  
% function phi = eval_ortho_matrix_param(n,p,A,B)
% parametrisation locale autour de B

global generalstiefel
if isempty(generalstiefel)
    generalstiefel=0;
end
if n==1 && ~generalstiefel && nargin<4
% HYPERSHEPE
    phi = zeros(1,p-1);
    phi(1)=acos(A(1));
    k=2;
    for k=2:p-2
        alpha = prod(sin(phi(1:k-1)));
        if alpha~=0
            phi(k)  =acos(A(k)/alpha);  
        else
            break
            warning('arret dans la recherche des parametres')
        end
    end
    if isempty(k)
        k=2;
    end

    if k==p-2 && A(end-1)~=0
        phi(end)=atan2(A(end),A(end-1));
    else
        phi(end)=0;
    end

    phi(end)=atan2(A(end),A(end-1));
    for k=p-2:-1:1
        phi(k)=atan2(A(k+1)/cos(phi(k+1)),A(k));
    end

else
    error('pas programme')
% VARIETE DE STIEFEL
    d = n*p - n*(n+1)/2;
    I = zeros(p,p);
    I(1:n*p)=1:n*p;
    I = tril(I,-1);
    I = nonzeros(I(:));
    S = zeros(p,p);
    S(I)=phi;
    S = S-S';
    A = expm(S);

    if nargin==4
        B = B';    
        if size(B,1)~=p
            error('la matrice locale n''a pas la bonne taille') 
        end
        if n<p
            [Q,R]=qr(B);
            B = [B,Q(:,n+1:p)];
            A = B*A*[eye(n);zeros(p-n,n)];
        else
            A = B*A;
        end
    else
        A = A*[eye(n);zeros(p-n,n)];    
    end
    A = A';

%elseif n==2 & p==2
%if phi<0    
%A = [cos(phi),sin(phi);sin(phi),-cos(phi)];
%else
%A = [cos(phi),sin(phi);-sin(phi),cos(phi)];    
%end

end


return

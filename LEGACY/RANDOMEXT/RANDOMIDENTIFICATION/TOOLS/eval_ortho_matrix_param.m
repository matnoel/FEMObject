function A = eval_ortho_matrix_param(n,p,phi,B)
%  A = eval_ortho_matrix_param(n,p,phi)
%  evaluation d'une matrice orthogonale de taille n-by-p parametree
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
% function A = eval_ortho_matrix_param(n,p,phi,B)
% parametrisation local autour de B

global generalstiefel
if isempty(generalstiefel)
    generalstiefel=0;
end
if n==1 & ~generalstiefel & nargin<4
% HYPERSHEPE
    A = zeros(n,p);

    A(1:p-1)=cos(phi);A(p)=1;
    for k=1:p-1
        A(k+1:end) = A(k+1:end)*sin(phi(k)); 
    end

else
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

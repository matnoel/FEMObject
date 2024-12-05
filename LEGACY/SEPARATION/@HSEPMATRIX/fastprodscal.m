function w = fastprodscal(u,v,A)
% Produit scalaire rapide :
%   marche pour (u,v) vecteurs,
%   fastprodscal(u,v)   = <u,v>   
%   fastprodscal(u,v,A) = <u,v>_A

if nargin==2     %   fastprodscal(u,v)   = <u,v>   
    ua=u.alpha;
    va=v.alpha;
    uF=u.F;
    vF=v.F;
    
%     dim=[1:u.dim];
%     w =zeros(u.m*v.m,1);
%     for i=1:u.m
%         for j=1:v.m
%             I=(i-1)*v.m+j;
%             w(I) = ua(i)*va(j);
%             for k=dim
%                 w(I) = w(I)*  fastprodscal( uF{i,k},vF{j,k} );
%             end
%         end
%     end
%     w=sum(w);
    
    [I,J]   = ind2sub([u.m,v.m],1:u.m*v.m);
    w       = cellfun(@fastprodscal,uF(I,:),vF(J,:));
    w       = sum(prod([w (ua(I).*va(J))'],2));
    
else           %   fastprodscal(u,v,A) = <u,v>_A
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % APPARAMENT LE PLUS RAPIDE :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w=fastprodscal(u,A*v);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Autre solution :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dim=[1:u.dim];
%     w=zeros(v.m*u.m,1);
%     ua=u.alpha;
%     va=v.alpha;
%     Aa=A.alpha;
%     for i=1:u.m
%         for j=1:v.m
%             I=(i-1)*v.m+j;
%             for ia=1:A.m
%                 wI = ua(i)*va(j)*Aa(ia);
%                 for k=dim
%                     wI = wI*  fastprodscal( u.F{i,k},v.F{j,k},A.F{ia,k} );
%                 end
%                 w(I)=w(I)+wI;
%             end
%         end
%     end
%     w=sum(w);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Autre solution :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dim=[1:u.dim];
%     [I,J,K]    = ind2sub([u.m,v.m,A.m],1:u.m*v.m*A.m);
%     w          = cellfun(@fastprodscal,u.F(I,dim),v.F(J,dim),A.F(K,dim));
%     w          = sum(prod([w (u.alpha(I).*v.alpha(J).*A.alpha(K))'],2));
end
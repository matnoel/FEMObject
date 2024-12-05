function R = RADIALMATRIX(V,L,sV,sL,D)
% function R = RADIALMATRIX(V,L,sV,sL)
% R= V1 tensorial L1 + V2 tensorial L2 + ... +  Vm tensorial Lm  (produit tensoriel)
% s : taille de V   s=[n1,n2]
% V : double (n1*n2)-by-m    ou  m cells contenant n1-by-n2 double
% L : double m-by-n3 ou m cells contenant (1-by-n3 double) ou autres objets
%                       se comportant comme des double (avec memes operations)
%
% function R = RADIALMATRIX(V,L,sV,sL,D)
% R= D(1,1)* V1 tensorial L1 + D(1,2)* V1 tensorial L2 + ...
% D : array of size (m,m)    identite par defaut

if nargin==1 & isa(V,'RADIALMATRIX')
    R = V ;    

elseif nargin>=2

    if nargin<3
        V = MULTIMATRIX(V,s);
    end

    if isa(L,'cell')
        L=vertcat(L{:});
    end

    R.m = size(V,2);
    V = MULTIMATRIX(V,s);

    if length(V)~=size(L,1)
        error('dimensions of V and L must match')      
    end

    R.V = V;
    R.L = L;

    if nargin==4
        R.D = D;
    else
        R.D=speye(R.m);
    end

    R= class(R,'RADIALMATRIX');
end


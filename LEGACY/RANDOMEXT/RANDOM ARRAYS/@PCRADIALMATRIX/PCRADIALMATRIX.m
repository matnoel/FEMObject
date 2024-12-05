function R = PCRADIALMATRIX(varargin)
% function R = PCRADIALMATRIX(V,sV,L)
% R= V1 tensorial L1 + V2 tensorial L2 + ... +  Vm tensorial Lm  (produit tensoriel)
% V : double (n1*n2)-by-m    ou  m cells contenant n1-by-n2 double ou
%  MULTIMATRIX de taille n1-by-n2
% L : PCMATRIX contenant les m variables aleatoires
%
% function R = PCRADIALMATRIX(V,sV,L,D)
% R= D(1,1)* V1 tensorial L1 + D(1,2)* V1 tensorial L2 + ...
% D : array of size (m,m)    identite par defaut
%
% function R = PCRADIALMATRIX(sV,PC)
% initialisation m=0
% sV : taille de la matrice
% PC : POLYCHAOS

if nargin==0
    R.V = MULTIMATRIX();
    R.L = PCMATRIX();
    R.m=[];
    R.DLmasse = cell(1,0);
    R.D=[];    

    R= class(R,'PCRADIALMATRIX',POLYCHAOS());
    superiorto('PCMATRIX')
elseif nargin==1 && isa(varargin{1},'PCRADIALMATRIX')
    R = varargin{1};

elseif nargin==2
    sV=varargin{1};
    PC = POLYCHAOS(varargin{2});
    R.V = MULTIMATRIX(sparse(prod(sV),0),sV,[0,1]);
    R.L = sparse([0,1],PC);
    R.m=0;
    R.DLmasse = cell(1,0);
    R.D=speye(0,0);    

    R= class(R,'PCRADIALMATRIX',PC);
    superiorto('PCMATRIX')
elseif nargin>=2
    V=varargin{1};
    if isa(V,'cell') || isa(V,'double')    
        sV=varargin{2};
        V=MULTIMATRIX(V,sV);
    elseif ~isa(V,'MULTIMATRIX')
        error('mauvais argument')
    end
    L = varargin{3};
    if isa(L,'cell')
        L = vertcat(L{:});
    end
    if isa(L,'PCRADIALMATRIX') 
        L = expand(L);    
    end
    if  ~isa(L,'PCMATRIX')
        error('L must be a PCMATRIX')
    end
    L=L(:);

    if length(V)~=size(L,1)
        error('number of L and V functions must equal')      
    end

    R.V = V;
    R.L = sparse(L);
    R.m = length(V);
    R.V = reshapem(R.V,[R.m,1]);
    PC = getPC(R.L);
    R.DLmasse=cell(1,0);
    if nargin==4
        R.D = varargin{4};
    else
        R.D=speye(R.m);
    end

    R= class(R,'PCRADIALMATRIX',PC);
    superiorto('PCMATRIX')

else
    error('bad arguments')

end


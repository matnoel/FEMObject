function R = TIMERADIALMATRIX(varargin)
% function R = TIMERADIALMATRIX(V,sV,L)
% R= V1 tensorial L1 + V2 tensorial L2 + ... +  Vm tensorial Lm  (produit tensoriel)
% V : double (n1*n2)-by-m    ou  m cells contenant n1-by-n2 double ou
%  MULTIMATRIX de taille n1-by-n2
% L : TIMEMATRIX contenant les m variables aleatoires
%
% function R = TIMERADIALMATRIX(V,sV,L,D)
% R= D(1,1)* V1 tensorial L1 + D(1,2)* V1 tensorial L2 + ...
% D : array of size (m,m)    identite par defaut
%
% function R = TIMERADIALMATRIX(sV,T)
% initialisation m=0
% sV : taille de la matrice
% T : TIMEMODEL

if nargin==0
    R.V = MULTIMATRIX();
    R.L = TIMEMATRIX();
    R.m=[];
    R.DLmasse = cell(1,0);
    R.D=[];    

    R= class(R,'TIMERADIALMATRIX',TIMEMODEL());
    superiorto('TIMEMATRIX')
elseif nargin==1 && isa(varargin{1},'PCRADIALMATRIX')
    R = varargin{1};
elseif nargin==2
    sV=varargin{1};
    T = gettimemodel(varargin{2});
    R.V = MULTIMATRIX(sparse(prod(sV),0),sV,[1,0]);
    R.L = zeros([0,1],T);
    R.m=0;
    R.DLmasse = cell(1,0);
    R.D=speye(0,0);    

    R= class(R,'TIMERADIALMATRIX',T);
    superiorto('TIMEMATRIX')
elseif nargin>=2
    V=varargin{1};
    if isa(V,'cell') || isa(V,'double')    
        sV=varargin{2};
        V=MULTIMATRIX(V,sV);

    elseif isa(V,'FEELEMFIELD')
        sV=[size(V,1),1];
    elseif ~isa(V,'MULTIMATRIX')
        error('mauvais argument')
    end
    L = varargin{3};
    if isa(L,'cell')
        L = vertcat(L{:});
    end
    if isa(L,'TIMERADIALMATRIX') 
        L = expand(L);    
    end
    if  ~isa(L,'TIMEMATRIX')
        error('L must be a TIMEMATRIX')
    end
    L=L(:);

    if ~isa(V,'FEELEMFIELD') && length(V)~=size(L,1)
        error('number of L and V functions must equal')      
    end

    R.V = V;
    R.L = sparse(L);
    R.m = length(V);
    T = gettimemodel(R.L);
    R.DLmasse=cell(1,0);
    if nargin==4
        R.D = varargin{4};
    else
        R.D=speye(R.m);
    end

    R= class(R,'TIMERADIALMATRIX',T);
    superiorto('TIMEMATRIX')

else
    error('bad arguments')

end


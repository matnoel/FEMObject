function R = PCRADIAL(V,L,D)
% function R = PCRADIAL(V,L)
% R= V1 tensorial L1 + V2 tensorial L2 + ... +  Vm tensorial Lm  (produit tensoriel)
% V : cell array of objects
% L : cell array of PCARRAY
%
% function R = PCRADIAL(V,L,D)
% R= D(1,1)* V1 tensorial L1 + D(1,2)* V1 tensorial L2 + ...
% D : array of size (m,m)    identite par defaut
if nargin==1 & isa(V,'POLYCHAOS')
    R.V=cell(1,0);
    R.L=cell(1,0);
    R.m=0;
    R.Lmasse = cell(1,0);
    R.D=speye(0,0);    

    R= class(R,'PCRADIAL',POLYCHAOS(V));

    superiorto('MYDOUBLE')
elseif nargin>=2

    if ~isa(V,'cell')
        V={V};
    end
    if ~isa(L,'cell')
        L = {L};
    end
    if ~isa(L{1},'PCARRAY')
        error('L must be a PCARRAY')
    end

    if length(V)~=length(L)
        error('number of cells in L and V must match')      
    end

    R.V = V;
    R.L = L;
    R.m = length(V);

    PC = getPC(R.L{1});


    R.Lmasse=cell(1,0);
    if nargin==3
        R.D = D;
    else
        R.D=speye(R.m);
    end

    R= class(R,'PCRADIAL',PC);
    superiorto('MYDOUBLE','PCARRAY')
end


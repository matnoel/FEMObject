function H = POLYFEND(nbdim,Z,varargin)
% H = POLYFEND(N)
% base elements finis L2 sur [0,1]^N, partition 
% polyn�mes orthonorm�s pour la mesure unitaire
% nbdim nombre de dimensions stochastiques
%
% H = POLYFEND(N,Z)
% Z : cellule de taille nbelem
% Z{j}.way : chemin de l'element
% Z{j}.order : ordre de l'element
%       si decoupage octree isotrope 


if nargin==1
    H = struct();
    H.e{1}.way = ones(1,nbdim);
    H.e{1}.order = 0;%zeros(0,nbdim);
    H.n = 1;
    h = cell(1,nbdim);
    h(:) = {POLYFE([0,1],[],'elem')};
    H.j = ones(1,nbdim);
    H.state = zeros(H.n,1);
    H.statepoint = []; 
    R = RANDPOLYS(h{:});
    H = class(H,'POLYFEND',R);

elseif isa(Z,'cell')
    H=struct();
    H.e = Z;
    H.n = length(Z);
    nbdim = size(Z{1}.way,2);
    H.j = zeros(H.n,nbdim);
    state = zeros(H.n,1);
    h = cell(1,nbdim);
    for k=1:nbdim
        A = zeros(H.n,2);
        for i=1:H.n
            A(i,:) = calc_xi(Z{i},k,nbdim);
        end
        [A,l,j] = unique(A,'rows');
        h{k} = POLYFE(A,[],'elem');
        H.j(:,k) = j; 
    end

    H.state = state;
    H.statepoint = [];
    R = RANDPOLYS(h{:});
    H = class(H,'POLYFEND',R);

end

function [xi] = calc_xi(E,n,nbdim)
% n : dimension selectionnee

if length(E.order)==1
    xi = sum((E.way(:,n)-1).*repmat((1/2).^(1:E.order)',1,1),1);
    xi = [xi;xi+(1/2)^(E.order)];
else

    xi = sum((E.way-1).*((1/2).^(E.order)),1);
    if size(E.order,1)==0
        xi = [xi;xi+(1/2).^zeros(1,nbdim)];    
    else
        xi = [xi;xi+(1/2).^(E.order(end,:))];
    end

    xi = xi(:,n);
end

return

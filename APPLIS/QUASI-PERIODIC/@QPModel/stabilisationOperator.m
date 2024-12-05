function op = stabilisationOperator(model,factor,connectivity,weightsOperator)
% op = stabilisationOperator(model,factor,connectivity,weightsOperator)

if nargin < 4
    weightsOperator = [] ;
end

if islogical(connectivity) % then it is a periodicity indication
    if isscalar(connectivity)
        connectivity = connectivity([1 1]);
    end
    connectivity = {mesoConnectivity(model,1,connectivity(1)) ; ...
        mesoConnectivity(model,2,connectivity(2))} ;
end

order = getOrder(model) ;
cellNum = getCellNum(model) ;
cModel = getCellModel(model) ;
cSz = getCellSize(model) ;
edges = getedges(getCellDomain(model)) ; % indexed as 4 1 2 3
op = [] ;
for d = 1:2
    % Basic operator
    conn = connectivity{d};
    [i,j] = find(conn) ;
    if isempty(i) ; continue ; end % safety (and saves time)
    cells = [formatIndex(order,cellNum,i) formatIndex(order,cellNum,j)] ;
    ccells = cell(4,1) ;
    % transpose connectivity for edge normal to direction d
    ccells{d} = cells(:,[order:2*(order-1) 1:order-1]) ;
    % do not transpose connectivity for edge normal to direction -d
    ccells{d+2} = cells ;
    op0 = bilinFormOperator(model,[0 0 0],ccells,factor) ;
    
    % Scale by face measure
    op0 = op0/cSz(mod(d,2)+1) ; % orthogonal direction
    
    % op0 is (in order 2 format)
    % {}^t\ki^d \otimes M_0^d
    %  + \ki^d \otimes M_0^{-d}
    % Note: face measure scaling is included in M_0^{±d} and M_1^{±d}
    
    % Apply stabilisation weights
    if ~isempty(weightsOperator)
        op0 = op0.*weightsOperator ;
    end % default stabilisation weight is 1
    
    % Lumping
    op1 = op0 ;
    for mu = 1:order-1
        sz = op1.space.sz(:,mu)' ;
        for c = 1:op1.space.dim(mu)
            op1.space.spaces{mu}{c} = spdiags(... % lump matrices
                sum(op1.space.spaces{mu}{c},1)',0,sz(1),sz(2)) ;
        end
    end
    
    % op1 is (in order 2 format)
    % l({}^t\ki^d \odot \omega) \otimes M_0^d
    %  + l(\ki^d \odot \omega) \otimes M_0^{-d}
    
    % Translation
    P = calc_P_edges(cModel,edges{mod(d,4)+1},edges{mod(d+2,4)+1}) ; % correct index
    op2 = op0 ;
    for c = 1:op2.space.dim(end)/2
        op2.space.spaces{end}{c} = P*op2.space.spaces{end}{c} ;
    end
    for c = 1+op2.space.dim(end)/2:op2.space.dim(end)
        op2.space.spaces{end}{c} = P'*op2.space.spaces{end}{c} ;
    end
    
    % op2 is (in order 2 format)
    % {}^t\ki^d \odot \omega \otimes M_1^d
    %  + \ki^d \odot \omega \otimes M_1^{-d}
    
    op = op + op1 - op2 ;
end
end
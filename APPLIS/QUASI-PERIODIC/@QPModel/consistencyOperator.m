function op = consistencyOperator(model,factor,connectivity,weightsOperator)
% op = consistencyOperator(model,factor,connectivity,weightsOperator)

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
edges = getedges(getCellDomain(model)) ; % indexed as 4 1 2 3
op = [] ;
for d = 1:2
    % Basic operator
    conn = connectivity{d} ;
    [i,j] = find(conn) ;
    if isempty(i) ; continue ; end % safety (and saves time)
    cells = [formatIndex(order,cellNum,i) formatIndex(order,cellNum,j)] ;
    ccells = cell(4,1) ;
    
    if isempty(weightsOperator)
        % transpose connectivity for edge normal to direction d
        ccells{d} = cells(:,[order:2*(order-1) 1:order-1]) ;
        % do not transpose connectivity for edge normal to direction -d
        ccells{d+2} = cells ;
        op0 = bilinFormOperator(model,[0 1 0],ccells,factor) ;
        op0 = .5*op0 ; % default average weight is 0.5
    else
        % For edge normal to direction d
        ccells{d} = cells ;
        op0 = bilinFormOperator(model,[0 1 0],ccells,factor) ;
        op0 = op0.*weightsOperator ;
        % Transpose mesoscopic operators
        for mu = 1:order-1
            op0.space.spaces{mu} = cellfun(@ctranspose,op0.space.spaces{mu}, ...
                'UniformOutput',false) ;
        end
        % For edge normal to direction -d
        ccells = cell(4,1) ;
        ccells{d+2} = cells ;
        op0 = op0 + bilinFormOperator(model,[0 1 0],ccells,factor).*weightsOperator ;
    end
    
    % op0 is (in order 2 format)
    % \sum_{n=1}^{r_K} {}^t (\ki^d \odot K^I_n \odot beta) \otimes N_0^d[K^Y_n]
    % + \ki^d \odot K^I_n \odot beta \otimes N_0^{-d}[K^Y_n]
    
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
    % \sum_{n=1}^{r_K} {}^t l(\ki^d \odot K^I_n \odot \beta) \otimes N_0^d[K^Y_n]
    % + l(\ki^d \odot K^I_n \odot \beta) \otimes N_0^{-d}[K^Y_n]
    
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
    % \sum_{n=1}^{r_K} {}^t \ki^q \odot K^I_n \odot \beta \otimes N_1^d[K^Y_n]
    % + l(\ki^q \odot K^I \odot \beta \otimes N_1^{-d}[K^Y_n]
    
    op = op + op1 - op2 ;
end
end
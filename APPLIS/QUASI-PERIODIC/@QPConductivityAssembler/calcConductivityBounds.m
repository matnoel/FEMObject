function conductivityBounds = calcConductivityBounds(assembler,method)
% conductivityBounds = calcConductivityBounds(assembler,method)
% Computes *square* of equivalency constants min and max from norm defined 
% by conductivity over model to L2 norm, over every ixcellModel for all i
% in cellNum.

if nargin == 1
    method = 1 ;
end

K = getConductivity(assembler) ;
cmodel = getCellModel(assembler) ;
cellNb = getCellNb(assembler) ;
distribution = getDistribution(assembler) ;
distribution = distribution(~cellfun(@isempty,distribution)) ; % safety
order = getOrder(assembler) ;
cellNum = getCellNum(assembler) ;

assert(prod(K.space.sz(1:order-1))==cellNb,'Conductivity incompatible with model')

switch method
    case 1
        % L^2 norm over the model
        m = BILINFORM(0,0,1) ;
        M = calc_matrix(m,cmodel) ;
        
        % Evaluate conductivity at one cell per phase
        numPhases = numel(distribution) ;
        cellIndices = cellfun(@(x) x(1),distribution(:)) ;
        cellIndices = formatIndex(order,cellNum,cellIndices) ;
        KEval = evalAtIndices(K,cellIndices,1:order-1) ;
        KEval = double(KEval)' ;
        
        Cmin = zeros(cellNb,1) ;
        Cmax = zeros(cellNb,1) ;
        for i = 1:numPhases
            khat = BILINFORM(0,0,KEval(:,i),0) ;
            Khat = calc_matrix(khat,cmodel) ;
            eigsPhase = eigs(Khat,M) ;
            eigsPhase = sort(eigsPhase) ; % advice from eigs documentation
            rep = distribution{i} ;
            Cmin(rep) = eigsPhase(1) ;
            Cmax(rep) = eigsPhase(end) ;
        end
    case 2 % not recommended
        C = 0 ;
        for i = 1:K.space.dim(end)
            max_i = max(K.space.spaces{end}(:,i)) ;
            C = max(C,max_i) ;
        end
        Cmax = C*ones(cellNb,1) ;
        Cmin = C*ones(cellNb,1) ;
end

conductivityBounds = [Cmin,Cmax] ;

if ~all(all(conductivityBounds>0))
    warning('Not all conductivityBounds are strictly positive.')
end
end

%TODO: chose definitely a method
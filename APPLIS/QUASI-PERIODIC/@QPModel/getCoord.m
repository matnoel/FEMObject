function x = getCoord(model,varargin)
% x = getCoord(model,varargin)

%TODO: keep best performing method (random test)
if ischarin('alt',varargin)
    x = getCoordAlt(model) ;
    return
end

compress = ~ischarin('noCompression',varargin) ;
if compress
    tolerance = getcharin('tolerance',varargin,getTolSVD(model)) ;
    maxRank = getcharin('maxRank',varargin,getCellNb(model)) ;
    Compressor = Truncator('tolerance',tolerance,'maxrank',maxRank) ;
else
    warning('getCoord: Compression is recommended') ;
end


y = getCellCoord(model) ;
L = getCellSize(model) ;
domainDim = getDim(model) ;
nbCell = getCellNb(model) ;
cellNum = getCellNum(model) ;
order = getOrder(model) ;
if order == 2 
    cellsSub = ind2NDsub(cellNum,(1:nbCell)') ;
end
excepted = getcharin('except',varargin,[]) ;
excepted = formatIndex(order,cellNum,excepted) ;

% Initialization
x = cell(1,domainDim) ;

for d = 1:domainDim
    x_space = cell(order,1) ;
    switch order
        case 2
            for i = 1:cellNum(d)
                x_space{1}(:,i) = double(cellsSub(:,d)==i) ;
                x_space{2}(:,i) = y(:,d)+(i-1)*L(d) ;
            end
            x_space{1}(excepted,:) = 0 ;
        case 3
            for i = 1:cellNum(d)
                for j = setdiff(1:domainDim,d)
                    x_space{j}(:,i) = ones(cellNum(j),1) ;
                    if ~isempty(excepted)
                        currentExcept = excepted(:,d)==i ;
                        currentExcept = excepted(currentExcept,j) ;
                        x_space{j}(currentExcept,i) = 0 ;
                    end
                end
                x_space{d}(:,i) = full(sparse(i,1,1,cellNum(d),1)) ;
                x_space{end}(:,i) = y(:,d) + (i-1)*L(d) ;
            end
    end
    x_space = TSpaceVectors(x_space) ;
    x_core = DiagonalTensor(ones(min(x_space.dim),1),order) ;
    % Could be optimized with no diagonal core
    x{d} = TuckerLikeTensor(x_core,x_space) ;
    
    % SVD
    if compress
        x{d} = Compressor.truncate(x{d}) ; % returns a CanonicalTensor
        x{d} = TuckerLikeTensor(x{d}.core,x{d}.space) ;
    end
end

end

function x = getCoordAlt(model)
% x = getCoordAlt(model)

% Build coordinates of full discontinuous domain
% (copied from calcTransferDiscontinuous)
cellNb = getCellNb(model) ;
cellNum= getCellNum(model);
cellCoord = getCellCoord(model) ;
coordDiscon = repmat(cellCoord,cellNb,1) ;
[i,j] = ind2sub(cellNum,(1:cellNb)');
shift = kron([i j]-1,ones(size(cellCoord,1),1)) ;
coordDiscon = coordDiscon + shift ;

% Extract continuous coordinates and exploit "tensorize" cell recursion
coordDiscon = mat2cell(coordDiscon,size(coordDiscon,1),ones(1,size(coordDiscon,2))) ;
x = tensorize(model,coordDiscon) ;
end
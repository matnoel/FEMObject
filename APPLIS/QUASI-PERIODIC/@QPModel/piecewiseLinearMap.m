function [map,mapGrad,detMapGrad] = piecewiseLinearMap(model,dilatation,limit)
% [map,mapGrad,detMapGrad] = piecewiseLinearMap(model,dilatation,limit)

% Dilatation format
cN = getCellNum(model) ;
if isnumeric(dilatation)
    dilatation={dilatation};
end
if numel(dilatation)==1
    dilatation=[dilatation {ones(cN(2)+1,1)}] ;
end
assert(all(cellfun(@numel,dilatation)==cN+1),...
    'Dilatation input has incorrect size')

% Limit format
if size(limit,1)==1
    limit = [limit ; -1 Inf] ;
end

% 
cellCoord = getCellCoord(model) ;
delta = max(limit,[],2)-min(limit,[],2) ;
mesoOffsetFun = @(d) delta(d)*((1:cN(d))'-1) + ...
    .5*(1-delta(d))*( cumsum([0;dilatation{d}(1:end-2)]) + ...
    cumsum([0;dilatation{d}(2:end-1)]) ) ;
mapFun = @(x,d) [ones(size(x,1),1)... % 1 everywhere
    min(x(:,d),limit(d,1)) ... % slope 1 on first band
    min(delta(d),max(0,x(:,d)-limit(d,1)))... % slope 1 on second band
    max(0,x(:,d)-limit(d,2))] ; % slope 1 on third band
mapGradFun = @(x,d) double( [x(:,d)<limit(d,1) ... % 1 on first band
    (x(:,d)>=limit(d,1) & x(:,d)<=limit(d,2)) ... % 1 on second band
    x(:,d)>limit(d,2)] ) ; % 1 on third band

% Loop on directions
map = [] ;
mapGrad = [] ;
detMapGrad = TuckerLikeTensor.ones(model.tensorSize) ;
for d = 1:2
    
    if all(dilatation{d}==1)
        % Map is coord in current direction and 0 in orthogonal one
        mapd = getCoord(model) ;
        mapd = mapd{d} ;
        mapd.space.spaces{end} = kron(mapd.space.spaces{end},sparse(d,1,1,2,1)) ;
        mapd = updateAllProperties(mapd) ;
        map = map + updateAllProperties(mapd) ;
        
        % Map gradient is 1 in current direction, 0 in other
        ts = model.tensorSize() ;
        ts(end) = 2*ts(end) ;
        mapGd = TuckerLikeTensor.ones(ts) ;
        mapGd.space.spaces{end}((3-d):2:end,:) = 0 ;
        mapGrad = mapGrad + mapGd ;
        continue
    end
    
    % Assemble temporary tensor for current direction
    mapMicro = kron(mapFun(cellCoord,d),sparse(d,1,1,2,1)) ;
    mapd = distributeMicro(model,...
        repmat({(1:prod(cN))'},1,size(mapMicro,2)),mapMicro) ;
    mapGMicro = kron(mapGradFun(cellCoord,d),sparse(d,1,1,2,1)) ;
    mapGd = distributeMicro(model,...
        repmat({(1:prod(cN))'},1,size(mapGMicro,2)),mapGMicro) ;
    
    % Apply dilatations
    mesoOffset = mesoOffsetFun(d) ;
    if getOrder(model)==3
        mapd.space.spaces{d}(:,1) = mesoOffset ;
        mapd.space.spaces{d}(:,2) = dilatation{d}(1:end-1) ;
        mapd.space.spaces{d}(:,4) = dilatation{d}(2:end) ;
        mapGd.space.spaces{d}(:,1) = dilatation{d}(1:end-1) ;
        mapGd.space.spaces{d}(:,3) = dilatation{d}(2:end) ;
    else
        [subs(:,1), subs(:,2)] = ind2sub(cN,(1:prod(cN))') ;
        for i = 1:cN(d) % loop on rows/columns
            subsdi = subs(:,d)==i ; % cells on row/column i
            mapd.space.spaces{1}(subsdi,1) = mesoOffset(i) ;
            mapd.space.spaces{1}(subsdi,2) = dilatation{d}(i) ;
            mapd.space.spaces{1}(subsdi,4) = dilatation{d}(i+1) ;
            mapGd.space.spaces{1}(subsdi,1) = dilatation{d}(i) ;
            mapGd.space.spaces{1}(subsdi,3) = dilatation{d}(i+1) ;
        end
    end
    
    % Sums
    map = map + mapd ;
    if nargout>2
        detMapGd = mapGd ;
        detMapGd.space.spaces{end} = detMapGd.space.spaces{end}(d:2:end,:) ;
        detMapGrad = detMapGrad.*updateAllProperties(detMapGd) ;
    end
    mapGrad = mapGrad + mapGd ;
end

mapGrad = toOperator(mapGrad) ;

end
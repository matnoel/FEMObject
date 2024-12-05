function op = bilinFormOperator(model,deriv,cells,factor)
% op = bilinFormOperator(model,deriv,cells,factor)

%% Preprocessing

if nargin < 4
    factor = [] ;
    if nargin < 3
        cells =  [] ;
        if nargin < 2
            deriv = [0 0 0] ;
        end
    end
end

if size(deriv,2) < 3
    deriv = deriv(:)' ;
    if numel(deriv) == 1
        deriv = repmat(deriv,1,2) ;
    end
    deriv = [deriv 0] ;
end

if isempty(cells)
    cells = repmat(formatIndex(getOrder(model),getCellNum(model),...
        (1:getCellNb(model))'),1,2) ;
end

if isempty(factor)
    factor = 1 ;
end

order = getOrder(model) ;
microModel = getCellModel(model) ;
factor = formatFactor(model,factor) ;
if iscell(factor)
    op = toOperator(factor{1}) ;
    factorRankMicro = factor{1}.space.dim(end) ;
elseif isa(factor.space,'TSpaceVectors')
    op = toOperator(factor) ;
    factorRankMicro = factor.space.dim(end) ;
elseif isa(factor.space,'TSpaceOperators')
    op = factor ;
    factorRankMicro = factor.space.dim(end) ;
end

if isnumeric(cells)
    %% Assembling assembling on (sub)domains (dim = d)
    
    % Mesoscopic
    if size(cells,2) == order-1
        cells = [cells cells] ;
    end
    mesoInd = mesoIndicatorT(model,cells) ;
    
    % Microscopic
    bf = setfree(BILINFORM(deriv(1),deriv(2),1,deriv(3)),0) ;
    for r = 1:factorRankMicro
        if iscell(factor)
            newk = cellfun(@(t)t.space.spaces{end}(:,min(r,t.space.dim(1)))...
                ,factor,'UniformOutput',false) ;
        elseif isa(factor.space,'TSpaceVectors')
            newk = factor.space.spaces{end}(:,r) ;
            if numel(newk)==2*getNbCellDoF(model) % vector field
                newk = {newk(1:2:end) ; newk(2:2:end)} ;
            end
        elseif isa(factor.space,'TSpaceOperators') % tensor field
            newk = factor.space.spaces{end}{r} ;
            ksz = size(newk) ;
            sub1 = (1:2:ksz(1))' ;
            sub2 = (2:2:ksz(1))' ;
            ind11 = sub2ind(ksz,sub1,sub1) ;
            ind12 = sub2ind(ksz,sub1,sub2) ;
            ind21 = sub2ind(ksz,sub2,sub1) ;
            ind22 = sub2ind(ksz,sub2,sub2) ;
            newk = {newk(ind11) newk(ind12) ; newk(ind21) newk(ind22)} ;
        end
        op.space.spaces{end}{r} = calc_matrix(setk(bf,newk),microModel) ;
    end
    
    op = mesoInd*op ; % does not commute
    
    % Previous code
    %     for o = 1:order-1
    %         opSpace{o} = cell(factor.space.dim(o),size(mesoInd,2));
    %         for i = 1:size(mesoInd,2)
    %             for j = 1:factor.space.dim(o)
    %                 opSpace{o}{i,j} = factor.space.spaces{o}(:,j).*mesoInd{o,i} ;
    %             end
    %         end
    %         opSpace{o} = opSpace{o}(:) ;
    %     end
    %     opSpace = TSpaceOperators(opSpace) ;
    %     % Core
    %     indRank = cellfun(@(x) size(x,2),mesoInd(1,:)) ;
    %     cumuRank = [0 cumsum(indRank)] ;
    %     opCore = zeros([cumuRank(end)*ones(1,order-1) 1]) ; % 1=opSpace.dim(end)/factor.space.dim(end)
    %     for d = find(indRank)
    %         subs = {(1+cumuRank(d)):cumuRank(d+1)} ;
    %         subs = [repmat(subs,1,order-1) {d}] ; % to use 2 or 3 (i.e. order) subscripts
    %         opCore(subs{:}) = eye(indRank(d)) ;
    %     end
    %     opCore = superkron(double(factor.core),opCore) ;
else
    %% Assembling on boundaries (dim = d-1)
    
    % Create boundary model, edges and normales
    boundary = create_boundary(microModel,'withparent') ;
    [edges,normales] = getedges(getCellDomain(model)) ;
    refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ;
    [~,orderInd] = ismember(refNormales,cell2mat(normales)','rows') ;
    normales = normales(orderInd) ;
    edges = edges(orderInd) ;
    
    if all(deriv(1:2)==0) % No derivative : BILINFORM
        bf = setfree(BILINFORM(0,0,1,deriv(3)),0) ;
        assblArgin = {} ;
    else % Derivative(s) : BILINFORMBOUNDARY
        bf = setfree(BILINFORMBOUNDARY(deriv(1),deriv(2),1,deriv(3)),0) ;
        assblArgin = {'parent',microModel} ;
    end
    notEmpty = find(~cellfun(@isempty,cells))' ;
    assert(~isempty(notEmpty),'No cell to assemble on.')
    edgeOp = op ;
    for c = notEmpty % process only non-empty
        % Mesoscopic
        if size(cells{c},2) == order-1
            cells{c} = [cells{c} cells{c}] ; % to get operators
        end
        mesoInd = mesoIndicatorT(model,cells{c}) ;
        
        % Microscopic
        cellEdgeModel = intersect(boundary,edges{c}) ;
        %         [~,edgeNumNodes] = intersect(microModel,edges{c}) ;
        for r = 1:factorRankMicro
            if iscell(factor)
                newk = cellfun(@(t)t.space.spaces{end}(:,min(r,t.space.dim(1))),factor,...
                    'UniformOutput',false) ;
                if deriv(1)~=deriv(2)
                    newk = {[newk{1,:}]*normales{c} ; ...
                        [newk{2,:}]*normales{c}} ;
                end
            elseif isa(factor.space,'TSpaceVectors')
                newk = factor.space.spaces{end}(:,r) ;
                if deriv(1)~=deriv(2)
                    newk = {normales{c}(1)*newk ; normales{c}(2)*newk } ;
                end
            elseif isa(factor.space,'TSpaceOperators')
                newk = factor.space.spaces{end}{r} ;
                ksz = size(newk) ;
                sub1x2 = repmat(1:2:ksz(1),2,1) ;
                sub2x2 = repmat(2:2:ksz(1),2,1) ;
                ind1 = sub2ind(ksz,(1:ksz(1))',sub1x2(:)) ;
                ind2 = sub2ind(ksz,(1:ksz(1))',sub2x2(:)) ;
                newk = {newk(ind1) ; newk(ind2)} ;
            end
            bf = setk(bf,newk) ;
            edgeOp.space.spaces{end}{r} = calc_matrix(bf,cellEdgeModel,...
                assblArgin{:}) ;
        end
        
        if c == notEmpty(1)
            op = mesoInd*edgeOp ; % does not commute
        else
            op = op + mesoInd*edgeOp ;
        end
    end
end
end

function factor = formatFactor(model,factor)

if isnumeric(factor)
    switch numel(factor)
        case 1
            factor = factor*TuckerLikeTensor.ones(tensorSize(model)) ;
        case getCellNb(model)
            factor = kron(factor,ones(getNbCellDoF(model),1)) ;
            factor = tensorize(model,factor) ;
            % TODO code more cost efficient
        case getNbCellDoF(model)
            factor = distributeMicro(model,{(1:getCellNb(model))'},factor) ;
    end
elseif iscell(factor)
    for i = 1:numel(factor)
        factor{i} = formatFactor(model,factor{i}) ;
    end
else
    assert(isa(factor,'TuckerLikeTensor'),'Input factor type unknown')
end

end
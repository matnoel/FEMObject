classdef QPModel
    
    properties( SetAccess = public, GetAccess = public )
        order % tensor order
        cellNum % number of cells along each direction
        cellSize % reference cell dimensions
        cellModel % reference cell model
        cellDomain
        elementSize % average element size along each direction
        tolSVD
        spatialTolerance % tolerance on coordinates and distances in spatial domain
        verbose
        discon2con
        con2discon
        traceConstantL2
    end
    
    methods( Access = public )
        %% Constructor methods
        function model = QPModel(varargin)
            p = ImprovedInputParser();
            addParameter(p,'order',2); % required
            addParameter(p,'cellNum',[2 2]); % required
            addParameter(p,'cellSize',[1 1]); % required
            addParameter(p,'elementSize',[.1 .1]); % required
            addParameter(p,'tolSVD',1e-6); % required
            addParameter(p,'spatialTolerance',1e-10); % required
            addParameter(p,'verbose',true); % required
            addParameter(p,'cellModel',[]);
            addParameter(p,'cellDomain',[]);
            addParameter(p,'discon2con',[]);
            addParameter(p,'con2discon',[]);
            addParameter(p,'traceConstantL2',[]);
            parse(p,varargin{:});
            model = passMatchedArgsToProperties(p,model);
            model = updateProperties(model) ; % Build missing attributes
        end
        
        function cellModel = buildCellModel(model)
            cDomain = getCellDomain(model) ;
            if isempty(cDomain)
                cDomain = buildCellDomain(model) ;
                model = setCellDomain(model,cDomain) ;
            end
            nElem = round(getCellSize(model)./getElementSize(model)) ;
            nElem = max(1,nElem) ; % for safety
            nElem = num2cell(nElem) ;
            cellModel = mesh(cDomain,nElem{:});
            cellModel = createddlnode(cellModel,DDL('u'));
        end
        
        function cellDomain = buildCellDomain(model)
            domainDim = getDim(model) ;
            lowerLeft = zeros(1,domainDim) ;
            upperRight = getCellSize(model) ;
            cellDomain = DOMAIN(domainDim,lowerLeft,upperRight) ;
        end
        
        
        function f = scalarField(model,generator)
            % f = scalarField(model,generator)
            % Generate continuous scalar field
            % Default generator is @rand
            if nargin < 2
                generator = @rand ;
            end
            f = generator(getNbDomainDoF(model),1) ;
            f = tensorize(model,f) ;
        end
        
        function f = boundaryScalarField(model,generator)
            % f = boundaryScalarField(model,generator)
            % Generate continuous scalar field whose support is boundary.
            % Default generator is @rand
            if nargin < 2
                generator = @rand ;
            end
            x = getDomainCoord(model) ;
            tol = 1e-15 ;
            loc = any(ismember(round(x/tol)*tol,[0 1]),2) ;
            f = zeros(getNbDomainDoF(model),1) ;
            f(loc) = generator(numel(find(loc)),1) ;
            f = tensorize(model,f) ;
        end
                
        %% "Get" methods
        
        % Primary
        
        function order = getOrder(model)
            order = model.order ;
        end
        
        function cellModel = getCellModel(model,varargin)
            cellModel = model.cellModel ;
            if ischarin('pbc',varargin)
               domain = getCellDomain(model) ;
               cellModel = addclperiodic(cellModel,getedge(domain,1),...
                   getedge(domain,3)) ;
               cellModel = addclperiodic(cellModel,getedge(domain,2),...
                   getedge(domain,4)) ;
            end
        end
        
        function cellNum = getCellNum(model)
            cellNum = model.cellNum ;
        end
        
        function cellSize = getCellSize(model)
            cellSize = model.cellSize ;
        end
        
        function cellDomain = getCellDomain(model)
            cellDomain = model.cellDomain ;
        end
        
        function elementSize = getElementSize(model)
            elementSize = model.elementSize ;
        end
        
        function tolSVD = getTolSVD(model)
            tolSVD = model.tolSVD ;
        end
        
        function spatialTolerance = getSpatialTolerance(model)
            spatialTolerance = model.spatialTolerance ;
        end
        
        function verbose = getVerbose(model)
            verbose = model.verbose ;
        end
        
        function traceConstantL2 = getTraceConstantL2(model)
            traceConstantL2 = model.traceConstantL2 ;
        end
        
        % Secondary
                
        function dim = getDim(model)
            dim = length(getCellSize(model)) ;
        end
        
        function nbDomainDoF = getNbDomainDoF(model)
            discon2Con = getDiscon2Con(model) ;
            nbDomainDoF = size(discon2Con,1) ;
        end
        
        function nbTotalDoF = getNbTotalDoF(model)
            nbTotalDoF = getCellNb(model)*getNbCellDoF(model) ;
        end
        
        function domainSize = getDomainSize(model)
            domainSize = getCellNum(model).*getCellSize(model) ;
        end
        
        function n = getCellNb(model)
            n = prod(getCellNum(model)) ;
        end
        
        % References
        
        function x = getCellCoord(model)
            x = getcoord(getnode(getCellModel(model))) ;
        end
        
        function n = getNbCellDoF(model)
            n = getnbddl(getCellModel(model)) ;
        end
        
        function discon2Con = getDiscon2Con(model)
            discon2Con = model.discon2con ;
        end
        
        function con2Discon = getCon2Discon(model)
            con2Discon = model.con2discon ;
        end
        
        %% "Set" methods
        
        function model = setOrder(model,newOrder)
            assert(ismember(newOrder,[2 3]),'Order can only be 2 or 3')
            model.order = newOrder ;
        end
        
        function model = setCellNum(model,newCellNum)
            model.cellNum = newCellNum ;
            model = updateProperties(model) ;
        end
        
        function model = setCellSize(model,newCellSize)
            model.cellSize = newCellSize ;
        end
        
        function model = setCellModel(model,newCellModel)
            if nargin == 1
                newCellModel = buildCellModel(model) ;
            end
            model.cellModel = newCellModel ;
            model = setTraceConstantL2(model) ;
            model = updateProperties(model) ; % constructor relies on this call
        end
        
        function model = setCellDomain(model,newCellDomain)
            if nargin == 1
                model = buildCellDomain(model) ;
            else
                model.cellDomain = newCellDomain ;
            end
        end
        
        function model = setElementSize(model,newElementSize)
            model.elementSize = newElementSize ;
        end
        
        function model = setTolSVD(model,newTolSVD)
            model.tolSVD = newTolSVD ;
        end        
        
        function model = setSpatialTolerance(model,newSpatialTolerance)
            model.spatialTolerance = newSpatialTolerance ;
        end
        
        function model = setVerbose(model,newVerbose)
            model.verbose = newVerbose ;
        end
        
        function model = setDiscon2Con(model,d2c)
            if nargin==1
                [d2c,c2d] = calcTransferDiscontinuous(model) ;
                model = setCon2Discon(model,c2d) ;
            end
            model.discon2con = d2c ;
        end
        
        function model = setCon2Discon(model,c2d)
            if nargin==1
                [d2c,c2d] = calcTransferDiscontinuous(model) ;
                model = setDiscon2Con(model,d2c) ;
            end
            model.con2discon = c2d ;
        end
        
        function model = setTraceConstantL2(model,newTraceConstantL2)
            if nargin == 1
                newTraceConstantL2 = calcTraceConstantL2(model) ;
            end
            model.traceConstantL2 = newTraceConstantL2 ;
        end
        
        %% External methods signatures
        
        corrected = applyCorrector(model,operator,corrector)
        
        v = basisVector(model,dir)
        
        op = bilinFormOperator(model,deriv,cells,factor)
        
        bCells = boundaryCells(model,cells,isExterior)
        
        [femModel,mesherTime] = buildFEModel(model,method)
        
        flux = calcBoundaryFlux(model,v,microOp,bCells,K)
        
        constant = calcTraceConstantL2(model)
        
        [discon2con,con2discon] = calcTransferDiscontinuous(model)
        
        [errorVal,diff] = compare2FE(model,tensQP,matFE,feCoord)
        
        c = component(model,t,dir)
        
        op = consistencyOperator(model,factor,connectivity,weightsOperator)
        
        t = createField(model,generator)
        
        tensor = distributeMicro(model,distribution,microVectors)

        x = extendTensor(model,x,xCells)
        
        x = fullCoord(model,i,y)
        
        x = getCoord(model,varargin)
        
        domainCoord = getDomainCoord(model)
        
        [] = ifprint(model,string)
        
        bounds = mesoBounds(model,t,dist)
        
        conn = mesoConnectivity(model,direction,isPeriodic)

        indicator = mesoIndicator(model,indices,isOperator)
        
        indicator = mesoIndicatorT(model,indices)
        
        [i,y] = mesoMicroCoord(model,x)
        
        penalty = penaltyBound(model,bounds,useSW,useAW,isPeriodic)
        
        [map,mapGrad,detMapGrad] = piecewiseLinearMap(model,dilatation,limit)
        
        [] = plot(model,f,varargin)
        
        [] = plotModes(model,tensor,gridSize)
         
        source = randomizeSource(model,maxSourceValue)
        
        x = restrictTensor(model,x,cells,orthogonalize)
        
        [op,srcOp] = sipOperator(model,K,source,isPeriodic,penalty)

        continuous = smooth(model,discontinuous)
        
        op = stabilisationOperator(model,factor,connectivity,weightsOperator)
        
        sModel = subModel(model,cells)
                
        discontinuous = sunder(model,continuous,coord)
        
        y = svdQP(model,m)
        
        [op,srcOp,penalty] = swipOperator(model,K,source,isPeriodic,penalty,useSW,useAW)

        v = switchBoundaryDoF(model,v,bCells,dir)
        
        t = tensorize(model,f,tolerance)
        
        tsz = tensorSize(model)
        
        f = untensorize(model,t,toSmooth)
        
        model = updateProperties(model)
        
        wOp = weightsOperator(model,bounds,type)
        
    end
    
    methods (Static)
        
        model = createRandom(varargin)
        
    end
end


classdef QPBC
    
    properties
        type
        loc
        value
        factor % either Robin factor or Nitsche penalty
    end
    
    methods( Access = public )
        
        function bc = QPBC(type,value,loc,factor)
            % bc = QPBC(type,value,loc,factor)
            if nargin < 4
                factor = [] ;
                if nargin < 3
                    loc = [] ;
                    if nargin < 2 || isempty(value)
                        value = [] ;
                    end
                end
            end
            bc.type = QPBC.formatType(type) ;
            if isempty(loc)
                tol = 1e-15 ;
                loc = @(x) any(ismember(round(x/tol)*tol,[0 1]),2) ;
            end
            bc.loc = QPBC.formatLocation(loc) ;
            if isempty(value)
                value = 0 ;
            end
            bc.value = QPBC.formatValue(value) ;
            assert(bc.type~=3||~isempty(factor),'Invalid Robin factor')
            bc.factor = factor ;
        end
        
        function t = eval(bc,model,asTensor)
            if nargin < 3
                asTensor = isa(model,'QPModel') ;
            end
            assert(numel(bc)==1,'QPBC arrays not supported')
            if isnumeric(model) && size(model,2)==2
                coord = model ;
                domainSz = max(coord) ;
            else
                coord = getDomainCoord(model) ;
                domainSz = getDomainSize(model) ;
            end
            coord = coord*diag(domainSz.^-1) ;
            t = bc.value(coord) ;
            iloc = bc.loc(coord) ;
            t(~iloc) = 0 ;
            if asTensor
                t = tensorize(model,t) ;
            end
        end
        
        function locCells = cells(bc,model)
            % locCells = cells(bc,model)
            
            % Ensure that this is a single bc
            assert(numel(bc)==1,'QPBC arrays not supported')
            
            % Get mesoscopic and microscopic coordinates.
            % The latter are mapped onto [0,1]^d
            coord = getDomainCoord(model) ;
            dSz = getDomainSize(model) ;
            x = coord*diag(dSz.^-1) ;
            cN = getCellNum(model) ;
            
            % Ensure that loc attribute matches only boundary cells
            locCells = unique(model.mesoMicroCoord(coord(bc.loc(x),:)),'rows') ;
            locCSub = formatIndex(3,cN,locCells) ;
            assert(all( ismember(locCSub(:,1),[1;cN(1)]) | ...
                ismember(locCSub(:,2),[1;cN(2)]) ),...
                'Invalid location: not all are boundary cells')
            
            % Sort cell subscripts in a 4-by-1 cell array 
            % (one edge per cell)
            locCells = {locCells(locCSub(:,1)==cN(1),:) ; ...
                locCells(locCSub(:,2)==cN(2),:) ; ...
                locCells(locCSub(:,1)==1,:) ; ...
                locCells(locCSub(:,2)==1,:) } ;
            
            % Test each cell.
            % This is to eliminate faces that are not involved in this BC
            % but were listed in locCells.
            % This may happen because a node (or at least its coordinates)
            % are shared across several cell faces.
            % E.g. on a 2x2 mesoscopic grid, considering a BC on edge 1
            % (i.e. of normale e_1), locCells would be so far
            % {[2;4] ; [4] ; [] ; [2]} 
            % whereas it should be
            % {[2;4]; [] ; [] ; []}.
            % The following block of code will remove those extra faces.
            spaTol = getSpatialTolerance(model) ;
            for e = find(~cellfun(@isempty,locCells(:)'))
                aOrthDir = mod(e,2)+1 ; % absolute orthogonal (to edge) direction
                locC = formatIndex(3,cN,locCells{e}) ;
                locC = locC(:,aOrthDir) ; % retrieve relevant cell subscripts
                toRemove = zeros(size(locC,1),1)==1 ;
                % List nodes that are on domain edge e
                edgeNodes = abs(x(:,setdiff([1 2],aOrthDir))...
                        -(1-sign(e-2.5))/2) < spaTol/cN(aOrthDir) ;
                % Test edge nodes on every cell face
                for i=1:size(locC,1)
                    % Restrict nodes that are interior to the cell face
                    interiorNodes = abs(x(:,aOrthDir)- ...
                        (locC(i)-.5)/cN(aOrthDir)) < (.5-spaTol)/cN(aOrthDir) ;
                    % N.B. interiorFaceNodes lists also nodes not on the
                    % domain edge (coordinates are considered along only
                    % one dimension). This list is only meant to be 
                    % intersected with edgeNodes.
                    %  Then, keep nodes that are on the domain edge and
                    % interior to the cell face.
                    testNodes = interiorNodes & edgeNodes ;
                    % Get test nodes' coordinates and test those against
                    % bc.loc to decide whether to keep this cell
                    toRemove(i) = ~any(bc.loc(x(testNodes,:))) ;
                end
                locCells{e}(toRemove,:) = [] ;
            end
            
            % Handle corner BC
            % If locCells is now empty, this must be because the BC applies
            % only on a corner.
            if all(cellfun(@isempty,locCells))
                % Identify domain corners in bc.loc and get their coordinates
                domainCorners = all(x==0 | x==1,2) ;
                locCoord = coord(bc.loc(x) & domainCorners,:) ;
                % Do previous cell identification again,
                locCells = unique(model.mesoMicroCoord(locCoord),'rows') ;
                locCSub = formatIndex(3,cN,locCells) ;
                locCells = {locCells(locCSub(:,1)==cN(1),:) ; ...
                    locCells(locCSub(:,2)==cN(2),:) ; ...
                    locCells(locCSub(:,1)==1,:) ; ...
                    locCells(locCSub(:,2)==1,:) } ;
            end
        end
        
        function [lhs,rhs] = operator(bc,model,K,applyValue,penalty)
            % [lhs,rhs] = operator(bc,model,K,applyValue,dirichletPenalty)
            % You should probably use operators instead.
            if nargin < 5
                penalty = [] ;
                if nargin < 4
                    applyValue = true ;
                end
            end
            assert(numel(bc)==1,'QPBC arrays not supported')
            
            if bc.type == 4
                lhs = [] ;
                rhs = [] ;
                return
            end
            
            locCells = bc.cells(model) ;
            locCells = cellfun(@(x)repmat(x,1,2),locCells,...
                'uniformoutput',false) ;
            
            switch bc.type
                case 1
                    if isempty(penalty)
                        if isempty(bc.factor)
                            penalty = 1 ;
                        else
                            penalty = bc.factor ;
                        end
                    end
                    M = penalty*bilinFormOperator(model,[0 0 0],locCells,1) ;
                    N = bilinFormOperator(model,[0 1 0],locCells,K) ;
                    lhs = M-N-N' ;
                    rhs = M-N' ;
                    
                case 2
                    rhs = bilinFormOperator(model,[0 0 0],locCells,1) ;
                    lhs = [] ;
                    
                case 3
                    lhs = bilinFormOperator(model,[0 0 0],locCells,bc.factor) ;
                    rhs = bilinFormOperator(model,[0 0 0],locCells,1) ;
            end
            
            X = getDomainCoord(model)*diag(getDomainSize(model).^-1) ;
            locInd = tensorize(model,double(bc.loc(X)),getTolSVD(model)/100) ;
            locInd = toQPOperatorIndicator(locInd,false) ; % right?
            if ~isempty(lhs)
                lhs = locInd.*lhs ; 
            end
            rhs = locInd.*rhs ;
            
            if applyValue
                val = bc.eval(model,true) ;
                if norm(val)<1e-14
                    rhs = [] ;
                else
                    rhs = rhs*val ;
                end
            end
        end
        
        function [lhs,rhs] = operators(bcs,model,K,applyValue,penalty)
            % [lhs,rhs] = operators(bcs,model,K,applyValue,penalty)
            % Assemble all LHS and RHS operators associated to provided
            % BCs. In case the solution is defined only up to a constant,
            % an operator \int_D uv is added to nullify it.
            if nargin < 5
                penalty = [] ;
                if nargin < 4
                    applyValue = true ;
                end
            end
            % Loop on BCs
            lhs=[];
            rhs=[];
            for n = 1:numel(bcs)
                [lhsk,rhsk] = bcs(n).operator(model,K,applyValue,penalty) ;
                lhs = lhs + lhsk ;
                rhs = rhs + rhsk ;
            end
            % In case of solution defined up to a constant. Apply penalty ?
            if undefinedConstant(bcs)
                lhs = lhs + bilinFormOperator(model,0) ;
            end
        end
        
        function bool = isPeriodic(bcs,dir)
            % bool = isPeriodic(bcs,dir)
            if nargin == 1
                dir = 1:2 ;
            end
            % Keep only periodic BC
            bcs = bcs([bcs(:).type]==4) ;
            % Build coordinates to be tested
            testX = (0:1e-3:1)' ;
            nX = numel(testX) ;
            testCoord = [] ;
            if any(dir==1)
                testCoord = [zeros(nX,1) testX ; ones(nX,1) testX] ;
            end
            if any(dir==2)
                testCoord = [testCoord ; ...
                    testX zeros(nX,1) ; testX ones(nX,1)] ;
            end
            % Test which coordinates have a PBC on it (no match test)
            flags = false ;
            for k = 1:numel(bcs)
                flags = flags | bcs(k).loc(testCoord) ;
            end
            % Return true iff all coordinates tested positive
            bool = all(flags) ;
        end
        
        function bool = undefinedConstant(bcs)
            % bool = undefinedConstant(bcs)
            % Returns true if provided BCs define solution only up to a
            % constant.
            types = [bcs(:).type] ;
            bool = all(ismember(types,[2 4])) ;
        end
    end
    
    methods (Static)
        
        function name = typeName(type)
            assert(isnumeric(type),'Type must be numeric');
            name = {'Dirichlet','Neumann','Robin','Periodic'} ;
            name = name{type} ;
        end
        
        function type = formatType(type)
            if ischar(type)
                if strcmpi(type,'dirichlet')
                    type = 1 ;
                elseif strcmpi(type,'neumann')
                    type = 2 ;
                elseif strcmpi(type,'robin')
                    type = 3 ;
                elseif strcmpi(type,'periodic') || strcmpi(type,'pbc')
                    type = 4 ;
                else
                    error('Unknow BC type')
                end
            elseif isscalar(type)
                assert(ismember(type,1:4),'QPBC type number must be 1, 2, 3 or 4')
            else
                error('Unknow BC type')
            end
        end
        
        function f = formatLocation(loc)
            if isa(loc,'function_handle')
                f = loc ;
            elseif isnumeric(loc) % edge numbers
                tol = 1e-15 ;
                f = @(~)false;
                if ismember(1,loc) % right edge
                    f = @(x)f(x) | 1-x(:,1)<tol ;
                end
                if ismember(2,loc) % top edge
                    f = @(x)f(x) | 1-x(:,2)<tol ;
                end
                if  ismember(3,loc) % left edge
                    f = @(x)f(x) | x(:,1)<tol ;
                end
                if ismember(4,loc) % bottom edge
                    f = @(x)f(x) | x(:,2)<tol ;
                end
            else
                error('Unknown BC location type')
            end
        end
        
        function f = formatValue(value)
            if isa(value,'function_handle')
                f = value ;
            elseif isnumeric(value)
                if isscalar(value)
                    f = @(x)value*ones(size(x,1),1) ;
                elseif any(size(value)==1)
                    f = @(x) value(:) ;
                else
                    error('Unknown BC value format')
                end
            else
                error('Unknown BC value type')
            end
        end
        
        function bc = random(type,generator)
            if nargin < 2
                generator = @rand ;
            end
            if nargin < 1 || isempty(type)
                type = randi(4) ;
            end
            value = [] ;
            loc = [] ;
            factor = [] ;
            if type~= 4
                value = @(x) generator(size(x,1),1) ;
                if type== 3
                    factor = generator(1,1) ;
                end
            end
            bc = QPBC(type,value,loc,factor) ;
        end
    end
end
classdef QPModelTest < matlab.unittest.TestCase
    
    properties
        maxIter = 100 ;
        testTol = 1e-10 ;
        tabCrit = @(a,b) max(max(abs(a-b))) ; % table error criterion
        matCrit = @(ref,rel) norm(rel-ref)/norm(ref) ; % matrix error criterion
    end
    
    methods (Static)
       
        function K = randomConductivity(model)
            % Hard-coded bounds
            maxValue = 100 ;
            minValue = 1 ;
            valGen = @() minValue-1+randi(maxValue-minValue+1) ;
            
            % Create patterns (random parameters)
            % uniform
            patterns(1) = struct('name','uniform','value',valGen(),...
                'size',[],'center',[],'offset',[]);
            % disc
            patterns(2) = struct('name','disc','value',valGen(),...
                'size',rand(1),'center',rand(1,2),'offset',[]) ;
            % rectangle'value',valGen(),
            patterns(3) = struct('name','rectangle','value',valGen(),...
                'size',rand(1,2),'center',rand(1,2),'offset',[]) ;
            % chevron
            patterns(4) = struct('name','chevron','value',valGen(),...
                'size',circshift([1 rand(1)],[0 randi(2)]),'center',rand(1,2),'offset',rand(1)) ;
            % bar
            patterns(5) = struct('name','bar','value',valGen(),...
                'size',circshift([1 rand(1)],[0 randi(2)]),'center',rand(1,2),'offset',[]) ;
            % cross
            patterns(6) = struct('name','cross','value',valGen(),...
                'size',rand(1,2),'center',rand(1,2),'offset',rand(1,2)) ;
            
            % Pattern table and distribution
            patternNb = numel(patterns) ;
            patternTable = rand(patternNb) ;
            dist = dealMultinomial(ones(1,patternNb)/patternNb,getCellNb(model));
            
            % Assemble
            cellCoord = getCellCoord(model) ;
            microK = drawCellPattern(cellCoord,patterns)*patternTable ;
            K = distributeMicro(model,dist,microK) ;
        end
        
    end
    
    methods (Test)
        
        function construction(testCase)
            noError = true ;
            for iter = 1:testCase.maxIter
                try
                    QPModel.createRandom('verbose',false) ;
                catch
                    noError = false ;
                    break
                end
            end
            testCase.verifyTrue(noError) ;
        end
        
        function fullCoord(testCase)
            model = QPModel.createRandom('verbose',false) ;
            cellSize = getCellSize(model) ;
            cellNum = getCellNum(model) ;
            cellNb = prod(cellNum) ;
            y = getCellCoord(model) ;
            [i,j] = ind2sub(cellNum,(1:cellNb)') ;
            mesoC = kron([i j],ones(size(y,1),1)) ;
            microC = repmat(y,cellNb,1) ;
            translation = (mesoC-1)*diag(cellSize) ;
            expX = microC + translation ;
            actX = fullCoord(model,mesoC,microC) ;
            testError = testCase.tabCrit(expX,actX) ;
            testCase.verifyLessThan(testError,testCase.testTol) ; 
        end
        
        function mesoMicroCoord(testCase)
            model = QPModel.createRandom('verbose',false) ;
            cellSize = getCellSize(model) ;
            cellNum = getCellNum(model) ;
            cellNb = prod(cellNum) ;
            y = getCellCoord(model) ;
            [i,j] = ind2sub(cellNum,(1:cellNb)') ;
            translation = ([i j]-1)*diag(cellSize) ;
            translation = kron(translation,ones(size(y,1),1)) ;
            expX = repmat(y,cellNb,1) + translation ;
            [acti,acty] = mesoMicroCoord(model,expX) ;
            actX = fullCoord(model,acti,acty) ;
            testError = testCase.tabCrit(expX,actX) ;
            testCase.verifyLessThan(testError,testCase.testTol) ; 
        end
        
        function transferDiscontinuous(testCase)
            model = QPModel.createRandom('verbose',false) ;
            cellCoord = getCellCoord(model) ;
            cellNum = getCellNum(model) ;
            cellNb = prod(cellNum) ;
            cellSize = getCellSize(model) ;
            [i,j] = ind2sub(cellNum,(1:cellNb)') ;
            translation = ([i j]-1)*diag(cellSize) ;
            translation = kron(translation,ones(size(cellCoord,1),1)) ;
            domainCoord = repmat(cellCoord,cellNb,1) + translation ;
            [~,~,expCon2Discon] = unique(domainCoord,'rows') ;
            actCon2Discon = getCon2Discon(model) ;
            con2DisconOK = all(actCon2Discon == expCon2Discon) ;
            expCon = rand(getNbDomainDoF(model),1) ;
            actCon = expCon(expCon2Discon) ;
            actCon = getDiscon2Con(model)*actCon ;
            discon2ConOK = all(actCon == expCon) ;
            testCase.verifyTrue(discon2ConOK && con2DisconOK) ;
        end
        
        function buildFEModel(testCase)
            model = QPModel.createRandom('verbose',false) ;
            feModel = buildFEModel(model) ;
            cellSize = getCellSize(model) ;
            cellNum = getCellNum(model) ;
            cellNb = prod(cellNum) ;
            y = getCellCoord(model) ;
            [i,j] = ind2sub(cellNum,(1:cellNb)') ;
            translation = ([i j]-1)*diag(cellSize) ;
            translation = kron(translation,ones(size(y,1),1)) ;
            expCoord = repmat(y,cellNb,1) + translation ;
            expCoord = expCoord(getDiscon2Con(model),:) ;
            actCoord = getcoord(getnode(feModel)) ;
            testError = testCase.tabCrit(sortrows(expCoord),...
                                         sortrows(actCoord)) ;
            testCase.verifyLessThan(testError,testCase.testTol) ;
        end
        
        function doubleVectorOrder2(testCase)
            model = QPModel.createRandom('verbose',false,'order',2) ;
            nbDoF = getNbDomainDoF(model) ;
            expVector = sunder(model,rand(nbDoF,1)) ;% duplicate faces DoF
            actVector = tensorize(model,expVector) ;
            actVector = doubleQP(actVector) ;
            vectorError = testCase.matCrit(expVector,actVector) ;
            testCase.verifyLessThanOrEqual(vectorError,testCase.testTol) ;
        end
        
        function doubleVectorOrder3(testCase)
            model = QPModel.createRandom('verbose',false,'order',3) ;
            nbDoF = getNbDomainDoF(model) ;
            expVector = sunder(model,rand(nbDoF,1)) ;% duplicate faces DoF
            actVector = tensorize(model,expVector) ;
            actVector = doubleQP(actVector) ;
            vectorError = testCase.matCrit(expVector,actVector) ;
            testCase.verifyLessThanOrEqual(vectorError,testCase.testTol) ;
        end
        
        function doubleOperatorOrder2(testCase)
            model = QPModel('order',2,'verbose',false,'cellNum',[2 3],...
                            'cellSize',[1.5 0.8],'elementSize',[0.1 0.2]) ; % curb DoF
            nbDoF = getNbDomainDoF(model) ;
            expOperator = sunder(model,rand(nbDoF)) ;% duplicate faces DoF
            actOperator = tensorize(model,expOperator) ;
            actOperator = doubleQP(actOperator) ;
            operatorError = testCase.matCrit(expOperator,actOperator) ;
            testCase.verifyLessThanOrEqual(operatorError,testCase.testTol) ;
        end
        
        function doubleOperatorOrder3(testCase)
            model = QPModel('order',2,'verbose',false,'cellNum',[2 3],...
                            'cellSize',[1.5 0.8],'elementSize',[0.1 0.2]) ; % curb DoF
            nbDoF = getNbDomainDoF(model) ;
            expOperator = sunder(model,rand(nbDoF)) ;% duplicate faces DoF
            actOperator = tensorize(model,expOperator) ;
            actOperator = doubleQP(actOperator) ;
            operatorError = testCase.matCrit(expOperator,actOperator) ;
            testCase.verifyLessThanOrEqual(operatorError,testCase.testTol) ;
        end
        
        function mesoBoundsOperator(testCase)
            model = QPModel.createRandom('verbose',false) ;
            K = QPModelTest.randomConductivity(model) ;
            Kop = K ;
            Kop.space.spaces{end} = kron(Kop.space.spaces{end},[1;1]) ;
            Kop.space = updateProperties(Kop.space) ;
            Kop = toOperator(Kop) ;
            ref = mesoBounds(model,K) ;
            rel = mesoBounds(model,Kop) ;
            err = testCase.tabCrit(ref,rel) ;
            testCase.verifyLessThanOrEqual(err,testCase.testTol) ;
        end
        
        function conductivityMatrix(testCase)
            model = QPModel.createRandom('verbose',false) ;
            K = QPModelTest.randomConductivity(model) ;
            [A1,b1] = swipOperator(model,K,true,[],true,true) ;
            Kop = K ;
            Kop.space.spaces{end} = kron(Kop.space.spaces{end},[1;1]) ;
            Kop.space = updateProperties(Kop.space) ;
            Kop = toOperator(Kop) ;            
            [A2,b2] = swipOperator(model,Kop,true,[],true,true) ;
            err = testCase.tabCrit(ref,rel) ;
            testCase.verifyLessThanOrEqual(err,testCase.testTol) ;
        end
        
%         function sipOperator(testCase)
%             pb =  QPDiffusionProblem.createRandom('verbose',true,'BC','PBC', ...
%                                                   'useAverageWeights',false, ...
%                                                   'useStabilisationWeights',false) ;
%             [ref,~,pb] = solve(pb,1:getOrder(pb)) ;
%             [A,b] = sipOperator(getModel(pb), getConductivity(pb), ...
%                                 getSource(pb), [true true]) ;
%             pb = updateGreedySolver(pb,A,b) ;
%             rel = solve(getGreedySolver(pb)) ;
%             err = testCase.matCrit(ref,rel) ;
%             testCase.verifyLessThan(err,testCase.testTol) ;
%         end
    end
    
end
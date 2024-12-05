classdef QPConductivityAssemblerTest < matlab.unittest.TestCase
    
    properties
        maxIter = 10
        testTol = 1e-10 ;
    end
    
    methods (Test)
        
        function construction(testCase)
            noError = true ;
            for iter = 1:testCase.maxIter
                try
                    QPConductivityAssembler.createRandom('verbose',false) ;
                catch
                    noError = false ;
                    break
                end
            end
            testCase.verifyTrue(noError) ;
        end
        
        function patterns(testCase)
            noError = true ;
            for iter = 1:testCase.maxIter
                try
                    QPConductivityAssembler.testPatterns('verbose',false) ;
                catch
                    noError = false ;
                    break
                end
            end
            testCase.verifyTrue(noError) ;
        end
        
        function assemble(testCase)
            assembler = QPConductivityAssembler.createRandom('verbose',false) ;
            assembler = assemble(assembler) ;
            actK = getConductivity(assembler) ;
            model = getModel(assembler) ;
            cellNb = getCellNb(model) ;
            f = getFields(assembler) ;
            dist = getDistribution(assembler) ;
            expK = zeros(cellNb*size(f,1),1) ;
            for p = 1:numel(dist)
                iphase = ismember((1:cellNb)',dist{p}) ;
                phase = kron(iphase,f(:,p)) ;
                expK = expK + phase ;
            end
            expK = smooth(model,expK) ;
            actK = untensorize(model,actK) ;
            tol = max(getTolSVD(assembler),getTolSVD(model)) ;
            Kerror = norm(actK-expK,'fro')/norm(expK,'fro') ;
            testCase.verifyLessThanOrEqual(Kerror,tol) ;
        end
        
        function conductivityBounds(testCase)
            assembler = QPConductivityAssembler.createRandom('verbose',false) ;
            order = getOrder(assembler) ;
            assembler = assemble(assembler) ;
            actBounds = getConductivityBounds(assembler) ;
            K = getConductivity(assembler) ;
            cmodel = getCellModel(assembler) ;
            cellNum = getCellNum(assembler) ;
            cellNb = getCellNb(assembler) ;
            m = BILINFORM(0,0,1) ;
            M = calc_matrix(m,cmodel) ;
            expBounds = zeros(cellNb,2) ;
            if order == 3
                [i,j] = ind2sub(cellNum,(1:cellNb)') ;
                cInd = [i,j] ;
            else
                cInd = (1:cellNb)' ;
            end
            K = double(evalAtIndices(K,cInd,1:order-1))' ;
            for i = 1:cellNb
                khat = BILINFORM(0,0,K(:,i),0) ;
                Khat = calc_matrix(khat,cmodel) ;
                eigsPhase = eigs(Khat,M) ;
                eigsPhase = sort(eigsPhase) ; % advice from eigs documentation
                expBounds(i,1) = eigsPhase(1) ;
                expBounds(i,2) = eigsPhase(end) ;
            end
            tol = getTolSVD(assembler) ;
            boundsError = max(max(abs((expBounds-actBounds)./expBounds))) ;
            testCase.verifyLessThan(boundsError,tol) ;
        end
        
    end
end
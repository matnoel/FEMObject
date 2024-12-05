function [ps,time] = assemble(ps,microOp,dGPenalty)
% [ps,time] = assemble(ps,microOp,dGPenalty)

if isempty(ps)
    warning('Empty patches. Will not assemble.')
    return
end

dbg = 4 ; % 1 for outside, 2 for inside, 3 for average and 4 for jump

clock = tic ;

K = getConductivity(ps) ;
if isempty(K)
    cA = assemble(getConductivityAssembler(ps)) ;
    ps = setConductivityAssembler(ps,cA) ;
    K = getConductivity(cA) ;
end
Kmax = getConductivityBounds(ps,'max') ;
model = getModel(ps) ;
patches = getPatches(ps) ;
pCells = getCells(ps) ;
if getUseCompression(ps)
    compressor = Truncator('tolerance',getTolSVD(model)) ;
end

rhsOp = cell(2,1) ;
for p = 1:numel(patches)
    switch getType(patches{p})
        case 1 % SWIP-like
            bcType = 4 ; % Cauchy
            isSIP = true ;
            useAW = true ;
            useSW = true ;
            penalty = dGPenalty ;
        case {2,3} % Neumann-Dirichlet, Lagrange
            bcType = 2 ; % Neumann
            isSIP = false ;
            useAW = false ;
            useSW = false ;
            penalty = dGPenalty ;
        otherwise
            error('Unexpected patch type.')
    end
    
    switch dbg
        case {1,2} % assemble on right side
            bcCells = boundaryCells(model,pCells{p}{1},dbg==1) ;
            bc = QPBC('type',bcType,'model',model,'cells',bcCells,...
                'isSIP',isSIP,'penalty',penalty,'useAverageWeights',useAW,...
                'useStabilisationWeights',useSW) ;
            [lhs,rhs] = assemble(bc,microOp,K,Kmax,0) ;
            for i = 2:numel(pCells{p})
                bcCells = boundaryCells(model,pCells{p}{i},dbg==1) ;
                bc = setCells(bc,bcCells) ;
                [newLHS,newRHS] = assemble(bc,microOp,K,Kmax,false) ;
                lhs = lhs + newLHS ;
                if bcType ~= 4
                    rhs = rhs + newRHS ;
                else
                    rhs{1} = rhs{1} + newRHS{1} ;
                    rhs{2} = rhs{2} + newRHS{2} ;
                end
            end
        case {3,4} % assemble on both side then average or substract
            for i = 1:numel(pCells{p})
                bcCellsOut = boundaryCells(model,pCells{p}{i},true) ;
                bcOut = QPBC('type',bcType,'model',model,'cells',bcCellsOut,...
                    'isSIP',isSIP,'penalty',penalty,'useAverageWeights',useAW,...
                    'useStabilisationWeights',useSW) ;
                [lhsOut,rhsOut] = assemble(bcOut,microOp,K,Kmax,0) ;
                bcCellsIn = boundaryCells(model,pCells{p}{1},false) ;
                bcIn = QPBC('type',bcType,'model',model,'cells',bcCellsIn,...
                    'isSIP',isSIP,'penalty',penalty,'useAverageWeights',useAW,...
                    'useStabilisationWeights',useSW) ;
                [lhsIn,rhsIn] = assemble(bcIn,microOp,K,Kmax,0) ;
                if dbg == 3 % average
                    newLHS = .5*lhsOut + .5*lhsIn ;
                    newRHS = .5*rhsOut + .5*rhsIn ;
                else % right - left, above - below
                    newLHS = lhsIn+lhsOut; % 0 anyway
                    rhsIn.space.spaces{1}{1} = -rhsIn.space.spaces{1}{1} ;
                    rhsIn.space.spaces{1}{2} = -rhsIn.space.spaces{1}{2} ;
                    rhsOut.space.spaces{1}{1} = -rhsOut.space.spaces{1}{1} ;
                    rhsOut.space.spaces{1}{2} = -rhsIn.space.spaces{1}{2} ;
                    newRHS = rhsOut + rhsIn ;
                end
                if i>1
                    lhs = lhs + newLHS ;
                    rhs = rhs + newRHS ;
                else
                    lhs = newLHS ;
                    rhs = newRHS ;
                end
            end
    end
    
    if getUseCompression(ps)
        if max(lhs.space.dim) > 1
            lhs = truncate(compressor,lhs) ;
        end
        if iscell(rhs)
            if max(rhs{1}.space.dim) > 1
                rhs{1} = truncate(compressor,rhs{1}) ;
            end
            if max(rhs{2}.space.dim) > 1
                rhs{2} = truncate(compressor,rhs{2}) ;
            end
        else
            if max(rhs.space.dim) > 1
                rhs = truncate(compressor,rhs) ;
            end
        end
    end
    if p == 1
        lhsOp = lhs ;
    else
        lhsOp = lhsOp + lhs ;
    end
    switch bcType
        case 2
            if p>1
                rhsOp{2} = rhsOp{2} + rhs ;
            else
                rhsOp{2} = rhs ;
            end
        case 4
            if p>1
                rhsOp{1} = rhsOp{1} + rhs{1} ;
                rhsOp{2} = rhsOp{2} + rhs{2} ;
            else
                rhsOp{1} = rhs{1} ;
                rhsOp{2} = rhs{2} ;
            end
    end
end
if isempty(rhsOp{1})
    rhsOp{1} = TuckerLikeTensor.zeros(tensorSize(model),tensorSize(model)) ;
end

ps = setLHSOperator(ps,lhsOp) ;
ps = setRHSOperator(ps,rhsOp) ;
time = toc(clock) ;
end
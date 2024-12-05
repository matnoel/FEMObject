function [sol,value,flux,output] = solve(ps,microOp,globalSource,globalSolution,tolerance,penalty)
% [sol,value,flux,output] = solve(ps,microOp,globalSource,globalSolution,tolerance,penalty)

patches = getPatches(ps) ;
pModel = getModel(ps) ;
pCells = getCells(ps) ;
output = cell(numel(patches),1) ;

% Loop on patches
for p = 1:numel(ps)
    % Preparations
    pNb = numel(pCells{p}) ; % number of patch instances
    source_p = cell(pNb,1) ;
    prevValue_p = cell(pNb,1) ;
    prevFlux_p = cell(pNb,1) ;
    
    % Pre-process QPPatch/solve input data
    for i = 1:pNb
        eBCells = boundaryCells(pModel,pCells{p}{i},1) ; % exterior boundary
        prevValue_p{i} = switchBoundaryDoF(pModel,globalSolution,eBCells) ;
        prevFlux_p{i} = calcBoundaryFlux(pModel,globalSolution,microOp,...
            eBCells,getConductivity(ps)) ;
        prevFlux_p{i} = switchBoundaryDoF(pModel,prevFlux_p{i},eBCells) ;
        subTargin = [mat2cell(pCells{p}{i},size(pCells{p}{i},1),...
            ones(size(pCells{p}{i},2),1)) {':'}] ;
        subTargin = cellfun(@unique,subTargin,'UniformOutput',false) ;
        source_p{i} = subTensor(globalSource,subTargin{:}) ;
        prevValue_p{i} = subTensor(prevValue_p{i},subTargin{:}) ;
        prevFlux_p{i} = subTensor(prevFlux_p{i},subTargin{:}) ;
    end
    
    % Resolution for all instances of current patch
    [sol_p,flux_p,output{p}] = solve(patches{p},source_p,prevValue_p,...
        prevFlux_p,tolerance,penalty) ;
    
    % Post-process QPPatch/solve output data
    for i = 1:pNb
        sol_p{i} = extendTensor(pModel,sol_p{i},pCells{p}{i}) ;
        flux_p{i} = extendTensor(pModel,flux_p{i},pCells{p}{i}) ;
        iBCells = boundaryCells(pModel,pCells{p}{i},0) ; % interior boundary
        if p==1 && i==1
            sol = sol_p{i} ;
            value = switchBoundaryDoF(pModel,sol_p{i},iBCells) ;
            flux = switchBoundaryDoF(pModel,flux_p{i},iBCells) ;
        else
            sol = sol + sol_p{i} ;
            value = value + switchBoundaryDoF(pModel,sol_p{i},iBCells) ;
            flux = flux + switchBoundaryDoF(pModel,flux_p{i},iBCells) ;
        end
    end
end

end
function [sol,flux,output] = solve(patch,source,value,flux,tolerance,penalty)
% [sol,flux,output] = solve(patch,source,value,flux,tolerance,penalty)

if nargin < 6
    penalty = 1 ;
end

if ~iscell(source)
    source = {source} ;
    value = {value} ;
    flux = {flux} ;
end
patchNb = numel(source) ;
assert(patchNb==numel(value)&&patchNb==numel(flux),'Bad input arguments.') ;

%% Build first problem

% For consistency between SWIP-like patch and diffusion assembler
useSW = true ;
useAW = true ;

% Create QP problem
pModel = getModel(patch) ;
pType = getType(patch) ;
cA = QPConductivityAssembler.fromConductivity(getConductivity(patch),pModel) ;
switch pType
    case 1 % SWIP-like
        bc = QPBC('type',4,'model',pModel,'value',{value{1};flux{1}},...
            'penalty',penalty,'isSIP',1,'useAverageWeights',useAW,...
            'useStabilisationWeights',useSW) ;
    case {2,3} % Dirichlet-Neumann, Lagrange
        bc = QPBC('type',1,'model',pModel,'value',value{1},...
            'penalty',1,'isSIP',false) ;
end
useMultiplier = pType == 3 ; % to be passed to solveFEM

dA = QPDiffusionAssembler('conductivityAssembler',cA,'BC',{bc},...
    'source',source{1},'useAverageWeights',useAW,...
    'useStabilisationWeights',useSW) ;
pb = QPDiffusionProblem('operatorAssembler',dA,'tolerance',tolerance) ;

%% Resolutions

sol = cell(patchNb,1) ;
flux = cell(patchNb,1) ;
switch getSolver(patch)
    case 1
        assert(pType~=3,'QP solver not implemented for Lagrange patch.')
        output = QPDiffusionProblem.qpOutputStructure(patchNb) ;
        [sol{1},output(1),pb] = solve(pb,1:getOrder(pb)) ;
        tic ;% Boundary fluxes (*outward*)
        bFlux = @(x) calcBoundaryFlux(pModel,x,getMicroOperators(pb),...
            boundaryCells(pModel),getConductivity(pb)) ;
        flux{1} = bFlux(sol{1}) ;
        output(1).time = output(1).time + toc ;
        % Loop on additional resolutions
        for n = 2:patchNb
            if pType==1
                bc = setValue(bc,value{n}) ;
            elseif pType==2
                bc = setValue(bc,{value{n};flux{n}}) ;
            end
            pb = setBC(pb,{bc}) ;
            pb = setSource(pb,source{n}) ;
            [sol{n},output(n),pb] = solve(pb,1:getOrder(pb)) ;
            tic ;
            flux{n} = bFlux(sol{n}) ;
            output(n).time = output(n).time + toc ;
        end
    case 2
        output = QPDiffusionProblem.feOutputStructure(patchNb) ;
        [sol{1},output(1)] = solveFEM(pb,useMultiplier) ;
        tic ;
        sol{1} = tensorize(pModel,sol{1}) ;
        % Boundary fluxes (*outward*)
        if useMultiplier
            flux{1} = tensorize(pModel,output(1).multiplier) ;
        else
            dA = assemble(getOperatorAssembler(pb)) ; % assemble micro operators
            bFlux = @(x) calcBoundaryFlux(pModel,x,getMicroOperators(dA),...
                boundaryCells(pModel),getConductivity(dA)) ;
            flux{1} = bFlux(sol{1}) ;
        end
        output(1).time(2) = output(1).time(2) + toc ; % add to assembling time
        % Loop on additional resolutions
        for n = 2:patchNb
            if pType==1
                bc = setValue(bc,value{n}) ;
            elseif pType==2
                bc = setValue(bc,{value{n};flux{n}}) ;
            end
            pb = setBC(pb,{bc}) ;
            pb = setSource(pb,source{n}) ;
            [sol{n},output(n)] = solveFEM(pb,useMultiplier) ;
            tic ;
            sol{n} = tensorize(pModel,sol{n}) ;
            if useMultiplier
                flux{n} = tensorize(pModel,output(n).multiplier) ;
            else
                flux{n} = bFlux(sol{n}) ;
            end
            output(n).time(2) = output(n).time(2) + toc ;
        end
end

end
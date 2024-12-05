function complete = assembleGlobalLocal(pb,glob,loc,fullModel)
% complete = assembleGlobalLocal(pb,glob,loc,fullModel)
% Combine a function over the global tensor product space and functions
% over local approximation spaces (one each) into a function over the
% complete space.
% _pb: QPDDProblem
% _glob: TuckerLikeTensor
% _loc: cell array of [double arrays or TuckerLikeTensors]
% _fullModel: MODEL
% _complete: TuckerLikeTensor if fullModel is empty, else double array

if nargin<4
    fullModel=[];
end

if ~iscell(loc)
    loc = {loc} ;
end

if isempty(fullModel) % Then global model is complete model. Process as tensors.
    
    % Ensure correct format
    if isnumeric(glob)
        if isscalar(glob)
            glob = glob*TuckerLikeTensor.ones(pb.model.tensorSize);
        else
            glob = tensorize(pb.model,glob) ;
        end
    end    
    % Set global solution to zero on patches
    nonPatchCells = setdiff(1:getCellNb(pb.model),pb.allPatchCells) ;
    globalRestricted = restrictTensor(pb.model,glob,nonPatchCells,...
        pb.compress) ;
    % Extend local solutions with zeros, as tensors
    localExtended = [] ;
    for i = 1:numel(loc)
        if isnumeric(loc{i})
            loc{i} = tensorize(pb.patch(i).model,loc{i}) ;
        end
        localExtended = localExtended + ...
            extendTensor(pb.model,loc{i},pb.patchCells{i}) ;
    end
    % Final assembling
    complete = globalRestricted + localExtended ;
    
else % Then fullModel is complete model (class MODEL). Process as double
    
    % Ensure correct format
    if isa(glob,'TuckerLikeTensor')
        globalRestricted = untensorize(pb.model,glob) ;
    elseif isnumeric(glob)
        if isscalar(glob)
            globalRestricted = glob(ones(getNbDomainDoF(pb.model),1)) ;
        else
            globalRestricted = smooth(pb.model,glob) ;
        end
    end
    % Get global solution outside patches
    patchLoc = untensorize(pb.model,mesoIndicatorT(pb.model,pb.allPatchCells))>0 ;
    globalRestricted = globalRestricted(~patchLoc) ;
    gCoord = getDiscon2Con(pb.model)*getDomainCoord(pb.model) ; % global continuous coord
    gCoord = gCoord(~patchLoc,:) ; % restrict to match globalRestricted
    % Reorder local solutions to fullModel node numbering (node assumed to
    % match inside patches), complete with zeros
    lCoord = cell(numel(loc),1) ;
    for i = 1:numel(loc)
        % Compute local coordinate offset with respect to global domain
        offset = min(pb.patchCells{i}) ; % index of lewer left patch cell
        offset = formatIndex(3,pb.model.cellNum,offset) ; % as subscripts
        offset = (offset-1)*diag(getCellSize(pb.model)) ; % then scale
        lCoord{i} = smooth(pb.patch(i).model,getDomainCoord(pb.patch(i).model)) ;
        lCoord{i} = [offset(1)+lCoord{i}(:,1) , offset(2)+lCoord{i}(:,2)];
        % Reorder (assuming match) and complete with zeros
        if isa(loc{i},'TuckerLikeTensor')
            loc{i} = untensorize(pb.patch(i).model,loc{i}) ;
        elseif isnumeric(loc{i})
            if isscalar(loc{i})
                loc{i} = loc{i}(ones(size(lCoord{i},1),1)) ;
            else
                loc{i} = smooth(pb.patch(i).model,loc{i}) ;
            end
        end
    end
    complete = cat(1,globalRestricted,loc{:}) ; % unsorted
    relCoord = cat(1,gCoord,lCoord{:}) ; % same order as "complete"
    refCoord = getcoord(getnode(fullModel)) ;
    assert(size(relCoord,1)==size(refCoord,1),'Mismatched dimensions.')
    complete = relativeSort(complete,relCoord,refCoord,pb.coordTol) ;
        
end
end
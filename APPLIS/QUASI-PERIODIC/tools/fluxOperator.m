function op = fluxOperator(model,loc,factor,nullOnBoundary,ptTol)
% op = fluxOperator(model,loc,factor,nullOnBoundary,ptTol)

if nargin < 5
    ptTol = 1e-10 ;
    if nargin < 4
        nullOnBoundary = true ;
        if nargin < 3
            factor = 1 ;
        end
    end
end

dofNb = getnbddl(model) ; % useful later
model = concatgroupelem(model,[],'enforce') ; % for calc_Rlocal below
% slow but assembling fails otherwise

% Create sub-model
% restricted to loc
locNodeNum = getnumber(model.node,find(loc)) ;
locModel = keepeleminnode(model,locNodeNum) ;
locModel = removeemptygroup(locModel);
numnodeelem = getnumnodeelem(locModel);
locModel = keepnode(locModel,numnodeelem);
locModel = calc_connec(locModel);

% Safety
if getnbelem(locModel)==0
    op = sparse(dofNb,dofNb) ;
    return
end

% Create boundary model of locModel
% but keep parent model relation
blModel = create_boundary(locModel,'withparent') ;
blModel = calc_connec(blModel) ;

% Get normale vector
% to each SEG2 element
normaleElem = calc_Rlocal(getgroupelem(blModel,1)) ;
% Store as element x dimension
normaleElem = double(permute(normaleElem,[3 2 1])) ;
normaleElem = normaleElem(:,[2 1]) ;
% Correct orientation for vertical elements
wrongSense = abs(normaleElem(:,2)) > 0.5 ; 
normaleElem(wrongSense,:) = -normaleElem(wrongSense,:) ;
% No, I don't understand why that works either.

% Convert to nodal field
% by averaging across connected elements
connecNodeElem = getfield(getconnec(blModel),'elem2node') ;
normaleNode = zeros(size(connecNodeElem,1),2) ;
blElemNum = getnumelem(blModel) ;
for i=1:size(normaleNode,1)
    % Get numeros of elements connected to node i
    connectedElemNum = find(connecNodeElem(i,:)~=0) ;
    iElem = ismember(blElemNum,connectedElemNum) ;
    % Average normales of those elements
    normaleNode(i,:) = mean(normaleElem(iElem,:),1) ;
    % Note: command below would likely work as well
    %     normaleNode(i,:) = mean(normaleElem(find(connecNodeEleme(i,:)~=0),:);1) ;
    % I.e. iElem is the same as connectedElemNum
end
% Normalisation
normaleNode = normaleNode./repmat(sqrt(sum(normaleNode.^2,2)),1,2) ;

% Extend to original model
% with zeros
normale = sparse(dofNb,2) ;
inode = getnumber(blModel.node) ;
normale(inode,:) = normaleNode ;

% Modify factor to include normale
factorAsCell = {factor(:,1).*normale(:,1) ; factor(:,end).*normale(:,2)} ;

% Assemble operator
op = calc_matrix(BILINFORMBOUNDARY(0,1,factorAsCell,0),blModel,'parent',model) ;
%
% For some reason (likely a bug), this assembles over model, with zeros on
% dof associated with nodes outside locModel. Consequently, the following
% commented lines are unnecessary (and return a "dimension mismatch" error).
%
% op = sparse(dofNb,dofNb) ;
% op(loc,loc) = calc_matrix(BILINFORMBOUNDARY(0,1,factorAsCell,0),blModel,...
%     'parent',locModel) ;

if nullOnBoundary
    relCoord = getcoord(getnode(model)) ;
    relCoord = relCoord-min(relCoord) ;
    relCoord = relCoord*diag(max(relCoord).^-1);
    boundaryLoc = abs(relCoord(:,1)) < ptTol | ...
        abs(relCoord(:,2)) < ptTol | ...
        abs(relCoord(:,1)-1) < ptTol | ...
        abs(relCoord(:,2)-1) < ptTol ;
    op(boundaryLoc,:) = 0 ;
end
end
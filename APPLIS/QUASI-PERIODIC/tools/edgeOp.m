function op = edgeOp(model,f)
% op = edgeOp(model,f)

if nargin < 2
    f = 1 ;
elseif iscell(f)
    f = horzcat(f{:}) ;
end

% Edges
coord = getcoord(getnode(model)) ;
domain = DOMAIN(getdim(model),min(coord),max(coord)) ;
[edges,normales] = getedges(domain) ;
refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ;
[~,orderInd]=ismember(refNormales,cell2mat(normales)','rows') ;
edges = edges(orderInd) ;
edgeNb = numel(edges) ;

% Build edges models
boundary = create_boundary(model) ;
edgeModel = cell(1,edgeNb) ;
for e=1:edgeNb
    edgeModel{e} = intersect(boundary,edges{e}) ;
end

edgeForm = BILINFORM(0,0,1,0) ;
edgeForm = setfree(edgeForm,0) ; % Unfree to apply no BC
op = cell(edgeNb,size(f,2)) ;
for r = 1:size(f,2) % For each pattern
    for e = 1:edgeNb % Per dim (in order) : N_0^ed N_0^(-ed)
        edgeForm = setk(edgeForm,f(:,r)) ;
        op{e,r} = calc_matrix(edgeForm,edgeModel{e}) ;
    end
end
end
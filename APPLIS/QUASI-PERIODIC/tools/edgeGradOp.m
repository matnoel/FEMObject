function op = edgeGradOp(model,f)
% op = edgeGradOp(model,f)

if nargin < 2
    f = 1 ;
elseif iscell(f)
    f = horzcat(f{:}) ;
end

% Normales and Edges
coord = getcoord(getnode(model)) ;
domain = DOMAIN(getdim(model),min(coord),max(coord)) ;
[edges,normales] = getedges(domain) ;
refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ;
[~,orderInd]=ismember(refNormales,cell2mat(normales)','rows') ;
normales = normales(orderInd) ;
edges = edges(orderInd) ;
edgeNb = numel(edges) ;

% Build edges models
boundary = create_boundary(model,'withparent') ;
edgeModel = cell(1,edgeNb) ;
for e=1:edgeNb
    edgeModel{e} = intersect(boundary,edges{e}) ;
end

edgeGradientForm = BILINFORMBOUNDARY(0,1,1,0) ;
edgeGradientForm = setfree(edgeGradientForm,0) ; % Unfree to apply no BC
op = cell(edgeNb,size(f,2)) ;
for r = 1:size(f,2) % For each pattern
    for e = 1:edgeNb % Per dim (in order) : N_0^ed N_0^(-ed)
        fn = {normales{e}(1)*f(:,r) ; normales{e}(2)*f(:,r) } ;
        edgeGradientForm = setk(edgeGradientForm,fn) ;
        op{e,r} = calc_matrix(edgeGradientForm,...
            edgeModel{e},'parent',model) ;
    end
end
end
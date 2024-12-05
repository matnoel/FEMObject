function traceConstant = traceConstant(cellModel)
% traceConstant = traceConstant(cellModel)

% Form over whole domain
M = calc_matrix(BILINFORM(1,1,1),cellModel) ;

% Edge models and form
coord = getcoord(getnode(cellModel)) ;
domain = DOMAIN(getdim(cellModel),min(coord),max(coord)) ;
edges = getedges(domain) ;
boundary = create_boundary(cellModel,'withparent') ;
edgeGradientForm = BILINFORMBOUNDARY(1,1,1,0) ;

traceConstant = 0 ;
for i = 1:4
    edgeModel = intersect(boundary,edges{i}) ;
    B = calc_matrix(edgeGradientForm,edgeModel,'parent',cellModel) ;
    eigenvalBM1 = eigs(B,M,1) ;
    if eigenvalBM1 < 0
        warning(['Strictly negative eigenvalue. Absolute value taken'])
    end
    traceConstant = max(traceConstant, sqrt(abs(eigenvalBM1))) ;
end
end
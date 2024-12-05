function constant = calcTraceConstantL2(model)
% C = calcTraceConstantL2(model)
% Merely calls function traceConstant on cell model

constant = traceConstant(getCellModel(model)) ;

end

% Legacy code
%
% version = 2 ;
% 
% cellModel = getCellModel(model) ;
% cellSize = getCellSize(model) ;
% cellMeasure = prod(cellSize) ;%sqrt(cellSize(1)^2 + cellSize(2)^2) ;
% 
% % L^2 norm on the domain
% m = BILINFORM(0,0,1) ;
% M = calc_matrix(m,cellModel) ;
% 
% switch version
%     case 1
%         % Create boundary model and set faces
%         cellDomain = getCellDomain(model) ;
%         Boundary = create_boundary(cellModel,'withparent') ;
%         [Edges,~] = getedges(cellDomain) ;
%         for e=1:4
%             cellModel = setfacet(cellModel,e,intersect(Boundary,Edges{e})) ;
%         end
%         traceConstant = 0 ;
%         for i=1:4
%             B = calc_matrix(m,getfacet(cellModel,i)) ;
%             newC = eigs(B,M,1) ;
%             newC = sqrt(newC*cellMeasure/getLength(Edges{i})) ;
%             traceConstant = traceConstant + newC ;
%         end
%         %TODO : is it possible to assemble over Boundary ?
%         
%     case 2
%         B = edgeOp(cellModel) ;
%         B = B{1}+B{2}+B{3}+B{4} ;
%         traceConstant = sqrt(eigs(B,M,1)*cellMeasure/min(cellSize)) ;
%         
%     case 3
%         M = calc_matrix(BILINFORM(1,1,1),cellModel) ;
%         B = edgeGradOp(cellModel) ;
%         traceConstant = 0 ;
%         for i = 1:numel(B)
%             traceConstant = max(traceConstant, sqrt(eigs(B{i},M,1))) ;
%         end
% end
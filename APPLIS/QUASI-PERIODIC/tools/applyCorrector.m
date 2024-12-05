function corrected = applyCorrector(operator,corrector,coord)
% corrected = applyCorrector(operator,corrector,coordinates)
% Apply corrector fields to get (apparent) homogenized quantity

corrected = zeros(2) ;
if isnumeric(operator)
    for i=1:2
        for j=1:2
            corrected(i,j) = (coord(:,i)+corrector(:,i))'*operator*coord(:,j) ;
        end
    end
    domainMeasure = prod(max(coord)) ;
elseif isa(operator,'AlgebraicTensor')
    for i=1:2
        for j=1:2
            corrected(i,j) = dot(coord{i}+corrector{i},operator*coord{j}) ;
        end
    end
    % Calculate domain measure
    uRi = coord{1}.sz(1:end-1) ; % upper right cell index/indices
    uRcoord = [double(evalAtIndices(coord{1},uRi,1:numel(uRi))) ; ...
        double(evalAtIndices(coord{2},uRi,1:numel(uRi)))] ;     
    domainMeasure = prod(max(uRcoord,[],2)) ;
else
    error('Not implemented')
end

% scaling
corrected = corrected/domainMeasure ;
end
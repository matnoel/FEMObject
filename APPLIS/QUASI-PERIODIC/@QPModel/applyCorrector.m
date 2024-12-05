function corrected = applyCorrector(model,operator,corrector)
% corrected = applyCorrector(model,operator,corrector)
% Apply corrector fields to get (apparent) homogenized quantity

coord = getCoord(model) ;
dim = getDim(model) ;

assert(dim==numel(corrector),'Inccorect number of correctors')

corrected = zeros(dim) ;
for i=1:dim
    for j=1:dim
        corrected(i,j) = dot(coord{i}+corrector{i},operator*coord{j}) ;
    end
end

% scaling
domainMeasure = prod(getDomainSize(model)) ;
corrected = corrected/domainMeasure ;
% TODO : why is scaling needed ?
end
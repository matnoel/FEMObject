function u= canonicalTensorOperator(u)
% function u = canonicalTensorOperator(u)
% conversion from SEPMATRIX to CanonicalTensor with TSpaceOperators

space = cell(size(u.F,2),1);
for i=1:size(u.F,2)
space{i} = u.F(:,i);
end
alpha = u.alpha;
u = CanonicalTensor(TSpaceOperators(space),alpha);
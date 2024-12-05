function u= canonicalTensor(u)
% function u = canonicalTensor(u)
% conversion from SEPMATRIX to CanonicalTensor with TSpaceVectors

space = cell(size(u.F,2),1);
for i=1:size(u.F,2)
space{i} = full(cell2mat(u.F(:,i)'));
end
alpha = u.alpha;
u = CanonicalTensor(TSpaceVectors(space),alpha);
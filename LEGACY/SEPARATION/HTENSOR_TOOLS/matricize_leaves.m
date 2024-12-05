function U = matricize_leaves(A)
% function U = matricize_leaves(A)

order = ndims(A);
U = cell(order,1);
r = rank(A);
r = r(A.dim2ind);
n = sqrt(size(A));
dim2ind = A.dim2ind;

for mu = 1:order
    U{mu} = cell(r(order),1);
    for i = 1:r(mu)
        U{mu}{i} = reshape(A.U{dim2ind(mu)}(:,i),n(mu),n(mu));
    end
end

end

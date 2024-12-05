function X = orth_basis(X)



tmp = cellfun(@(a)a(:) , X.A ,'UniformOutput',0);
tmp = [tmp{:}];

[L,D,U] = svd(tmp,'econ');

for k=1:size(L,2)
    X.A{k} = reshape(L(:,k),X.s);
end

X.Phi = D*U'*X.Phi;


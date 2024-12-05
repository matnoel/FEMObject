function X = compress(X,tol)

X = orth_basis(X);

if nargin<2
    tol = 1e-10;
end

% Cross approx of Phi : new function Phi, and evaluation point PivotElement
Phi = X.Phi;
[PivotValue,PivotElement,Rows,Cols,ifail] = CompleteACA(Phi',tol);
PivotElement=PivotElement(:,1);
n_eval = length(PivotElement);
Phi = (Cols/Cols(PivotElement,:))';

A=cell(n_eval,1);
parfor i=1:n_eval
    A{i} = eval_sparse(X,PivotElement(i)) ;
end



X = affine_matrix(A,Phi);


end
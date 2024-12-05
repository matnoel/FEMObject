function Z = mtimes_truncation(X,Y)

tol = 1e-10;

% Compute the Phi function of the product X*Y
[I,J]=ind2sub([Y.r,X.r],1:(X.r*Y.r));
Phi = Y.Phi(I,:).*X.Phi(J,:);

% Cross approx of Phi : new function Phi, and evaluation point PivotElement
[PivotValue,PivotElement,Rows,Cols,ifail] = CompleteACA(Phi',tol);
PivotElement=PivotElement(:,1);
n_eval = length(PivotElement);
Phi = (Cols/Cols(PivotElement,:))';

% Evaluation of the product on PivotElement
Ye = eval(Y,PivotElement);
Z_eval = cellfun( @(x) x*Ye , X.A,'UniformOutput',0 );
XPhi = X.Phi(:,PivotElement)';
XPhi = XPhi(:);
XPhi = sparse([1:length(XPhi)]',repmat([1:n_eval]',X.r,1),XPhi,length(XPhi), n_eval );

% Z_eval = [Z_eval{:}]*XPhi;
% Z_eval = mat2cell(Z_eval,size(Z_eval,1),ones(size(Z_eval,2),1))';

Z_eval=cell(n_eval,1);
parfor i=1:n_eval
    Xe = eval_sparse(X,PivotElement(i));
    Ye = eval_sparse(Y,PivotElement(i));
    Z_eval{i} = Xe*Ye;
end

Z = affine_matrix(Z_eval,Phi);

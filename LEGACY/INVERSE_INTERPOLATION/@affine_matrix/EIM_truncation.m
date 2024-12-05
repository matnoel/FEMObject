function X = EIM_truncation(X)
% Truncation of the affine matrix X based on the Casenave's method.
% +> EIM on X.Phi
% +> Evaluation of X on a (small) collection of point
% => Reduction of the rank X.r

tol = 1e-10;

Phi = X.Phi;

[PivotValue,PivotElement,Rows,Cols,ifail] = CompleteACA(Phi',tol);
PivotElement=PivotElement(:,1);
Phi = (Cols/Cols(PivotElement,:))';

X_eval=cell(length(PivotElement),1);
for i=1:length(PivotElement)
    X_eval{i} = eval(X,PivotElement(i));
end

X = affine_matrix(X_eval,Phi);


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

if n_eval>10
    
     Z_eval = cell(n_eval,1);
    parfor i=1:n_eval
        
        Yi = eval_sparse(Y,PivotElement(i));
        Z_eval{i} = eval(X*Yi,PivotElement(i));
        
    end
    
%     Y_eval = cell(n_eval,1);
%     parfor i=1:n_eval
%         Y_eval{i} = eval(Y,PivotElement(i));
%     end
%     Y_eval=[Y_eval{:}];
%     Z_eval = cell(X.r,1);
%     parfor i=1:X.r
%         Z_eval{i}=X.U{i}\(X.L{i}\(Y_eval* diag(X.Phi(i,PivotElement)) ));
%     end
%     
%     
%     s=size(Z_eval{1});
%     Z_eval = cellfun(@(z)z(:),Z_eval,'UniformOutput',0);
%     Z_eval = [Z_eval{:}];
%     Z_eval = sum(Z_eval,2);
%     Z_eval = reshape(Z_eval,s);
%     Z_eval = mat2cell(Z_eval,size(Z_eval,1),ones(1,size(Z_eval,2)));
    
    
    
else
    Z_eval = cellfun(@(ii)  eval(Y,ii) , num2cell(PivotElement),'UniformOutput',0);
    Z_eval = cellfun(@(Z_eval)  X*Z_eval , Z_eval,'UniformOutput',0);
    Z_eval = cellfun(@(Z_eval,ii) eval( Z_eval ,ii), Z_eval, num2cell(PivotElement),'UniformOutput',0);
end


%
%
% Ye = eval(Y,PivotElement);
% Z_eval = cellfun( @(x) x*Ye , X.A,'UniformOutput',0 );
% XPhi = X.Phi(:,PivotElement)';
% XPhi = XPhi(:);
% XPhi = sparse([1:length(XPhi)]',repmat([1:n_eval]',X.r,1),XPhi,length(XPhi), n_eval );
%
% Z_eval = [Z_eval{:}]*XPhi;
% Z_eval = mat2cell(Z_eval,size(Z_eval,1),ones(size(Z_eval,2),1))';

% Z_eval=cell(n_eval,1);
% for i=1:n_eval
%     Ye = eval(Y,PivotElement(i));
%     M = cellfun( @(x) x*Ye , X.A,'UniformOutput',0 );
%     M = cellfun( @(x) x(:) , M,'UniformOutput',0);
%     M = [M{:}]*X.Phi(:,PivotElement(i));
%     Z_eval{i} = M;
% end

Z = affine_matrix(Z_eval,Phi);

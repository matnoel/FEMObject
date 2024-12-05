function [U,S,V,err] = svdtruncate(X,tol,p,varargin)
% function [U,S,V,err] = svdtruncate(X,tol,p)
% function [U,S,V] = svdtruncate(X,tol,p)
% Performs a truncated economy-size singular value decomposition (SVD) of the matrix X
% such that Y = U*S*V' approximates X with relative precision tol
% in schatten-p norm (1<=p<=Inf, p=2 for Frobenius), p=2 by default
% if tol>=1, tol corresponds to the maximal rank of the SVD
% 
% function [s,err] = svdtruncate(X,tol,p)
% function s = svdtruncate(X,tol,p)
% Returns the first singular values of matrix X in descending order

if nargin==1 || isempty(tol)
    tol = [];
end
if nargin<3
    p = 2;
end

[U,S,V] = svd(full(X),'econ');
s = diag(S);

% if p==Inf
%     err = s/max(s);
% else
%     if verLessThan('matlab','8.2') % compatibility (<R2013b)
%         err = (flipdim(cumsum(flipdim(s.^p,1)),1)/sum(s.^p)).^(1/p);
%     else
%         err = (flip(cumsum(flip(s.^p)))/sum(s.^p)).^(1/p);
%     end
% end
% err = [err(2:end);0];
if p==Inf
    err = s/max(s);
    err = [err(2:end);0];
else
    err = (1-cumsum(s.^p)/sum(s.^p)).^(1/p);
end
if ~isempty(tol)
    if tol<1
        m = find(err<tol);
        if isempty(m)
            m = min(size(X));
        else
            m = min(m);
        end
    else
        m = min(tol,min(size(X)));
    end
else
    m = min(size(X));
end

U = U(:,1:m);
S = S(1:m,1:m);
V = V(:,1:m);
err = err(1:m);

if nargout<=2
    U = s(1:m);
    S = err(1:m);
end


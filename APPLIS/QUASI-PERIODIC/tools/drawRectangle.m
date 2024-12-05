function rectanglePattern = drawRectangle(X,size,relativeCenter)

if nargin < 3 || isempty(relativeCenter)
    relativeCenter = 0.5*[1 1];
end

if numel(size) == 1
    size = size*[1 1] ;
end

tol = 1e-6 ; % hard coded tolerance on coordinates

domainSize = max(X,[],1)-min(X,[],1) ;
C = domainSize.*relativeCenter ;

rectanglePattern = abs(X(:,1)-C(1)) - size(1)/2 < tol & ...
    abs(X(:,2)-C(2)) - size(2)/2 < tol;

end
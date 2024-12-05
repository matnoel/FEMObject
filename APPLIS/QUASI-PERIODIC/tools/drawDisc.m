function discPattern = drawDisc(X,radius,relativeCenter)

if nargin<3
    relativeCenter = [0.5 0.5];
end

tol = 1e-9 ; % hard coded tolerance on coordinates

domainSize = max(X,[],1)-min(X,[],1) ;
center = domainSize.*relativeCenter ;

discPattern = (X(:,1)-center(1)).^2 + (X(:,2)-center(2)).^2 - radius^2 < tol ;

end
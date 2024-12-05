function losangePattern = drawLosange(X,radii,relativeCenter)

if nargin<3
    relativeCenter = [0.5 0.5];
end

if isscalar(radii)
    radii = radii([1 1]);
end

tol = 1e-9 ; % hard coded tolerance on coordinates

domainSize = max(X,[],1)-min(X,[],1) ;
center = domainSize.*relativeCenter ;

losangePattern = abs(X(:,1)-center(1))/radii(1) + ...
    abs(X(:,2)-center(2))/radii(2) - 1 < tol ;

end
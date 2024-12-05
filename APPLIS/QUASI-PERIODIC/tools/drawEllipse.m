function ellipsePattern = drawEllipse(X,radii,relativeCenter,rotation)

if nargin<4
    rotation = 0 ;
    if nargin<3
        relativeCenter = [0.5 0.5];
    end
end

if isscalar(radii)
    radii = radii([1 1]) ;
end

tol = 1e-9 ; % hard coded tolerance on coordinates

domainSize = max(X,[],1)-min(X,[],1) ;
center = domainSize.*relativeCenter ;

% Translation
X = [X(:,1)-center(1) X(:,2)-center(2)] ;

% Rotation
R = [cos(rotation) -sin(rotation) ; sin(rotation) cos(rotation)] ;
X = X*R ;

ellipsePattern = (X(:,1)/radii(1)).^2 + (X(:,2)/radii(2)).^2 - 1 < tol ;

end
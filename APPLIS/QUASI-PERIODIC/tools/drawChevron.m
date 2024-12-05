function chevronPattern = drawChevron(X,width,offset,relativeCenter)

if nargin < 4 || isempty(relativeCenter)
    relativeCenter = .5*[1 1];
end

tol = 1e-9 ; % hard coded tolerance on coordinates

domainSize = max(X,[],1)-min(X,[],1) ;
C = domainSize.*relativeCenter ;
direction = find(ismember(width,domainSize),1) ;
assert(~isempty(direction),'Could not determine direction.')
orthDir = 1+mod(direction+2,2) ; % 2 if direction=1, 1 if direction=2
width = width(orthDir) ;

if nargin < 3
    % Set max offset with chevron within domain
    offset = 1 - (C(direction)+width)/domainSize(direction) ;
elseif isempty(offset)
    offset = 0 ;
end

chevron = offset*(1-abs(X(:,direction)/C(direction)-1)) ;
chevronPattern = (abs(X(:,orthDir) - (C(orthDir)+chevron) )) - width/2 < tol ;

end
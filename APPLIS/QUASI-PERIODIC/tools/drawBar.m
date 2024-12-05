function barPattern = drawBar(X,width,relativeCenter)
% function barPattern = drawBar(X,width,direction,relativeCenter)

% Draw as chevron with no offset
barPattern = drawChevron(X,width,0,relativeCenter) ;

% if nargin<4
%     relativeCenter = 0.5;
% end
% 
% domainSize = max(X,[],1)-min(X,[],1) ;
% orthDir = 1+mod(direction+2,2) ; % 2 if direction=1, 1 if direction=2
% center = domainSize(orthDir)*relativeCenter ;
% barPattern = abs(X(:,orthDir)-center) <= width/2 ;

end
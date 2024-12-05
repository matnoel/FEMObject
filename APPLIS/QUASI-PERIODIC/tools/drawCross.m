function crossPattern = drawCross(X,width,relativeCenter,offset)

if nargin < 4 || isempty(offset)
    offset = [0 0] ;
    if nargin < 3
        relativeCenter = 0.5*[1 1] ;
    end
end

% Draw as two chevrons
chevron1 = drawChevron(X,[width(1) 1],offset(1),relativeCenter) ;
chevron2 = drawChevron(X,[1 width(2)],offset(2),relativeCenter) ;

crossPattern = chevron1 | chevron2 ;

end
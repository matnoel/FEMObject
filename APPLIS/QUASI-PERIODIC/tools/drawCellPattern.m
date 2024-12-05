function patternField = drawCellPattern(coord,pattern)
% patternField = drawCellPattern(coord,pattern)

patternNb = numel(pattern) ;

if patternNb > 1 % then recurse
    patternField = zeros(size(coord,1),patternNb) ;
    for p = 1:patternNb
        patternField(:,p) = drawCellPattern(coord,pattern(p)) ;
    end
else
    switch pattern.name
        case 'uniform'
            patternField = ones(size(coord,1),1) ;
        case 'disc'
            patternField = drawDisc(coord,pattern.size,pattern.center) ;
        case 'rectangle'
            patternField = drawRectangle(coord,pattern.size,pattern.center) ;
        case 'chevron'
            patternField = drawChevron(coord,pattern.size,pattern.offset,pattern.center) ;
        case 'bar'
            patternField = drawBar(coord,pattern.size,pattern.center) ;
        case 'cross' 
            patternField = drawCross(coord,pattern.size,pattern.center,pattern.offset) ;
        case 'losange'
            patternField = drawLosange(coord,pattern.size,pattern.center) ;
        case 'ellipse'
            patternField = drawEllipse(coord,pattern.size,pattern.center,pattern.offset) ;
        otherwise
            error('Unknwon pattern name')
    end
    patternField = pattern.value*double(patternField) ;
end
end
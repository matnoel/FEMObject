function discontinuous = sunder(model,continuous,coord)
% discontinuous = sunder(model,continuous,coord)

if nargin == 3
    continuous = relativeSort(continuous,getDomainCoord(model),coord) ;
end

% Safety
if length(continuous) == getNbTotalDoF(model)
    discontinuous = continuous ;
    return
end

c2d = getCon2Discon(model) ;

szc = size(continuous) ;
if szc(1) > szc(2) 
    discontinuous = continuous(c2d,:) ;
elseif szc(1) == szc(2) % operator
    discontinuous = continuous(c2d,c2d) ;
else
    discontinuous = continuous(:,c2d) ;
end
end
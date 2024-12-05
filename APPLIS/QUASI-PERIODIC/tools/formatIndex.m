function I = formatIndex(order,cellNum,I)
% I = formatIndex(order,cellNum,I)

% Cell recursion
if iscell(I)
    for n = 1:numel(I)
        I{n} = formatIndex(order,cellNum,I{n}) ;
    end
    return
end

% Safety
if isempty(I) ; return ; end

% Ensure input format
sz = size(I) ;
if all(sz>2)
    error('Input cannot have more than two rows and more than two columns')
elseif sz(2)>2
    I = I' ;
    sz = size(I) ;
end

% Format according to order
if order == 2 && sz(2) == 2
    I = sub2ind(cellNum,I(:,1),I(:,2)) ;
elseif order == 3 && sz(2) == 1
    [i,j] = ind2sub(cellNum,I) ;
    I = [i j] ;
end

end
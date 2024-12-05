function x = restrictTensor(x,cells,cellNum,orthogonalise)
% x = restrictTensor(x,cells,cellNum,orthogonalise)
% Return restriction of x to cells in QP format sense : output has same
% domain but zero values outside cells.
%
% x : TuckerLikeTensor of order 2 or 3 (TSpaceVectors or TSpaceOperators)
% cells : column of cell indices or subscripts to which x will be
% restricted. If cellNum is not provided, must be formatted to suit order.
% cellNum ; [optional] domain mesoscopic grid size [n1 n2]. Provide to 
% format cells.
% orthogonalise : [optional] Boolean default false. Set to true to
% orthogonalise.

if nargin < 4
    orthogonalise = false ;
    if nargin < 3
        cellNum = [] ;
    end
end

if isempty(cells)
    x = [] ;
    return
end

if ~isempty(cellNum)
    cells = formatIndex(x.order,cellNum,cells) ;
end

if isa(x.space,'TSpaceOperators')
    for o = 1:x.order-1
        setZero = setdiff((1:x.space.sz(1,o))',cells(:,o)) ;
        for c = 1:x.space.dim(o)
            x.space.spaces{o}{c}(setZero,:) = 0 ;
            x.space.spaces{o}{c}(:,setZero) = 0 ;
        end
    end
    if orthogonalise
        warning('Orthogonalisation not implemented for TSpaceOperators.')
    end
else % is a TSpaceVectors
    for o = 1:x.order-1
        setZero = setdiff((1:x.space.sz(o))',cells(:,o)) ;
        x.space.spaces{o}(setZero,:) = 0 ;
    end
    if orthogonalise
        x = orth(x) ;
    end
end

% if orthogonalise
%     x = orth(x) ;
% end
%TODO: debug TSpace/orth for TSpaceOperators
x = updateAllProperties(x) ;

end


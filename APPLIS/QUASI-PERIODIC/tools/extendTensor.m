function x = extendTensor(x,cellNum,xCells)
% x = extendTensor(x,cellNum,xCells)
% Resizes x for a larger domain by extending with null values. x values are
% located at xCells.
%
% x : TuckerLikeTensor of order 2 or 3 (TSpaceVectors or TSpaceOperators)
% cellNum : Mesoscopic grid size [n1 n2] of new domain.
% xCells : Location of input x values in new domain.

% Ensure new size format matches order
mesoSize = formatIndex(x.order,cellNum,cellNum) ;

if isa(x.space,'TSpaceOperators')
    for o = 1:x.order-1
        new = zeros(mesoSize(o)) ;
        subs = unique(xCells(:,o)) ;
        for c = 1:x.space.dim(o)
            new(subs,subs) = x.space.spaces{o}{c} ;
            x.space.spaces{o}{c} = new ;
        end
    end
else % is a TSpaceVectors
    for o = 1:x.order-1
        subs = unique(xCells(:,o)) ;
        new = zeros(mesoSize(o),x.space.dim(o)) ;
        new(subs,:) = x.space.spaces{o} ;
        x.space.spaces{o} = new ;
    end
end

x = updateAllProperties(x) ;

end


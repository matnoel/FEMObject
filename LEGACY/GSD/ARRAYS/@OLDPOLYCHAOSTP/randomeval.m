function Hs = randomeval(PC,xi,RV)
% function Hs = randomeval(PC,xi,RV)

if nargin == 2
    RV = RANDVARS(PC);
    if (isa(xi,'cell') && length(xi)~=length(RV)) || (isa(xi,'double') && size(xi,2)~=length(RV))
        error(['on attend ' num2str(length(RV)) ' jeux de valeurs' ])
    end
    if isa(xi,'cell')
        xi = [xi{:}];
    end
end

if nargin==3
    xi = transfer(RANDVARS(RV),RANDVARS(PC),xi);
end

Hs = polyval(PC,xi);

function DAU = calc_rigitangu(S,q,U,varargin)
% DAU = calc_rigitangu(S,q,U,varargin)
% S : MODEL
% q : solution courante

q = unfreevector(S,q);
U = unfreevector(S,U);

DAU = calc_vector(S,@rigitangu,q,U,varargin{:});

if ~ischarin('nofree',varargin)
    DAU = freevector(S,DAU);
end

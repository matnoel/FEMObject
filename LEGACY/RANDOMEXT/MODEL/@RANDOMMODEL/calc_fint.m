function f = calc_fint(S,q,varargin)
% function f = calc_fint(S,q,varargin)
% S : MODEL
% q : solution courante

q = unfreevector(S,q);
f = calc_vector(S,@fint,q,varargin{:});

if ~ischarin('nofree',varargin)
    f = freevector(S,f);
end

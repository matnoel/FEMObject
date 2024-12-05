function K = calc_rigitang(S,q,varargin)
% function K = calc_rigitang(S,q,varargin)
% S : MODEL
% q : solution courante

q = unfreevector(S,q);
K = calc_matrix(S,@rigitang,q,varargin{:});

if ~ischarin('nofree',varargin)
    K = freematrix(S.BCOND,K);
end

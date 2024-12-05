function [K,b] = calc_rigi(S,varargin)
% function [K,b] = calc_rigi(S,varargin)
% calcul de la matrice de rigidite
% S : MODEL
% K : matrice de rigidite reduite
% b : dans le cas de conditions aux limites non-homogenes
%     si u=[u1;u2] avec u1 impose : b = K21*u1
%
% dans le cas ou S est aleatoire, appel de calc_rigipc(S,varargin{:})
% See also MODEL/calc_rigipc

if israndom(S)
    K = calc_rigipc(S,varargin{:});
else
    
    K = calc_matrix(S,@rigi,varargin{:});
    
    if nargout==2
        b = calc_nonhomogeneous_vector(S,K);
    end
    
    if ~ischarin('nofree',varargin)
        K = freematrix(S,K);
    end
    
end

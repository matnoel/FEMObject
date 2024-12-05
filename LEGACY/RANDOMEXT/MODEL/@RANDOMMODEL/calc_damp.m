function [C,b] = calc_damp(S,varargin)
% function [C,b] = calc_damp(S,varargin)
% calcul de la matrice d'amortissement
% S : MODEL
% C : matrice d'amortissement reduite
% b : dans le cas de conditions aux limites non-homogenes
%     si u=[u1;u2] avec u1 impose : b = C21*u1
%
% dans le cas ou S est aleatoire, appel de calc_damppc(S,varargin{:})
% See also MODEL/calc_damppc

if israndom(S)
    C = calc_damppc(S,varargin{:});
else
    
    if getnblevelsets(S)==0
        C = calc_matrix(S,@damp,varargin{:});
    else
        if ~isalevelset(S.ls)
            S = lseval(S);
        end
        C = calc_multimatrix(S,@mdampls,S.ls,varargin{:});
    end
    
    if ~ischarin('nofree',varargin)
        C = freematrix(S,C);
    end
    
end

function [M,b] = calc_mass(S,varargin)
% function [M,b] = calc_mass(S,varargin)
% calcul de la matrice de masse
% S : MODEL
% M : matrice de masse reduite
% b : dans le cas de conditions aux limites non-homogenes
%     si u=[u1;u2] avec u1 impose : b = M21*u1
%
% dans le cas ou S est aleatoire, appel de calc_masspc(S,varargin{:})
% See also MODEL/calc_masspc

if israndom(S)
    M = calc_masspc(S,varargin{:});
else
    
    if getnblevelsets(S)==0
        M = calc_matrix(S,@mass,varargin{:});
    else
        if ~isalevelset(S.ls)
            S = lseval(S);
        end
        M = calc_multimatrix(S,@massls,getlevelsets(S),varargin{:});
    end
    
    if nargout==2
        b = calc_nonhomogeneous_vector(S,M);
    end
    
    if ~ischarin('nofree',varargin)
        M = freematrix(S,M);
    end
    
end

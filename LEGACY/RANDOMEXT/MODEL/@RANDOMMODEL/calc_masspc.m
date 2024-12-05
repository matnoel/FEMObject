function [M,b] = calc_masspc(S,PC,varargin)
% function [M,b] = calc_masspc(S,PC,varargin)
% S : MODEL
% PC : POLYCHAOS
% M : matrice de masse reduite
% b : dans le cas de conditions aux limites non-homogenes
%     si u=[u1;u2] avec u1 impose : b = M21*u1

PC = getPC(PC);

if ~israndom(S)
    M = calc_mass(S).*one(PC);
else
    if getnblevelsets(S)==0
        M = calc_pcmatrix(S,PC,@masspc,varargin{:});
    else
        M = calc_pcmatrix(S,PC,@masslspc,getlevelsets(S),varargin{:});
    end
end

if nargout==2
    b = calc_nonhomogeneous_vector(S,M);
end

if ~ischarin('nofree',varargin)
    M = freematrix(S,M);
end

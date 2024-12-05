function [K,b] = calc_rigipc(S,PC,varargin)
% function [K,b] = calc_rigipc(S,PC)
% S : MODEL
% PC : POLYCHAOS
% K : matrice de rigidite reduite
% b : dans le cas de conditions aux limites non-homogenes
%     si u=[u1;u2] avec u1 impose : b = K21*u1

PC = getPC(PC);

if ~israndom(S)
    K = calc_rigi(S).*one(PC);
else
    if getnblevelsets(S)==0
        K = calc_pcmatrix(S,PC,@rigipc,varargin{:});
    else
        K = calc_pcmatrix(S,PC,@rigilspc,getlevelsets(S),varargin{:});
    end
end

if nargout==2
    b = calc_nonhomogeneous_vector(S,K);
end

if ~ischarin('nofree',varargin)
    K = freematrix(S,K);
end

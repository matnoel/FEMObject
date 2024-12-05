function [C,b] = calc_damppc(S,PC,varargin)
% function [C,b] = calc_damppc(S,PC,varargin)
% S : MODEL
% PC : POLYCHAOS
% C : matrice d'amortissement reduite
% b : dans le cas de conditions aux limites non-homogenes
%     si u=[u1;u2] avec u1 impose : b = C21*u1

PC = getPC(PC);

if ~israndom(S)
    C = calc_damp(S).*one(PC);
else
    if getnblevelsets(S)==0
        C = calc_pcmatrix(S,PC,@damppc,varargin{:});
    else
        C = calc_pcmatrix(S,PC,@damplspc,getlevelsets(S),varargin{:});
    end
end

if nargout==2
    b = calc_nonhomogeneous_vector(S,C);
end

if ~ischarin('nofree',varargin)
    C = freematrix(S,C);
end

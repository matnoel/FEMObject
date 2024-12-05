function K = calc_rigitangexpect(S,u,a,b,varargin)
% function K = calc_rigitangexpect(S,a,b,varargin)
% calcul de E(K*a*b) avec K matrice de rigidite
% S : MODEL
% a, b : PCMATRIX
%

if ~israndom(S) && ~israndom(u)
    K = calc_matrix(S,@rigitang,u,varargin{:});
    K = expect(K,a,b);
else
K = calc_multimatrix(S,@rigitangexpect,u,a,b,varargin{:});
end

K = freematrix(S.BCOND,K);


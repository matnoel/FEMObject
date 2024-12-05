function K = calc_rigiexpect(S,a,b,varargin)
% function K = calc_rigi(S,varargin)
% calcul de la matrice de rigidite
% S : MODEL
%

if ~israndom(S)
    K = calc_matrix(S,@rigi,varargin{:});
    K = expect(K,a,b);
else

K = calc_multimatrix(S,@rigiexpect,a,b,varargin{:});

end

K = freematrix(S.BCOND,K);


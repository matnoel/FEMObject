function K = calc_rigitangpc(S,q,PC,varargin)
% function K = calc_rigitangpc(S,q,PC,varargin)
% S : MODEL
% q : solution courante
% PC : POLYCHAOS

PC = getPC(PC);

q = unfreevector(S,q);
K = calc_pcmatrix(S,PC,@rigitangpc,q,varargin{:});
K = freematrix(S.BCOND,K);

if ~israndom(K)
    K = K.*one(PC);
end

function K = calc_rigitangpc(S,q,PC,varargin)
% function K = calc_rigitangpc(S,q,PC,varargin)
% S : MODEL
% q : solution courante
% PC : POLYCHAOS

PC = getPC(PC);

q = unfreevector(S,q);
K = calc_pcmatrix(S,PC,@rigitangpc,q,varargin{:});

if ~ischarin('nofree',varargin)
    K = freematrix(S.BCOND,K);
end

if ~israndom(K)
    K = K.*one(PC);
end

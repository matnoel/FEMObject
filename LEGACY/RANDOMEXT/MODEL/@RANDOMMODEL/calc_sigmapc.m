function s=calc_sigmapc(S,q,PC,varargin)

% function s=calc_sigmapc(S,q,PC)
% S : MODEL
% q : deplacement
% PC : POLYCHAOS
q=unfreevector(S,q);
if ischarin('node',varargin)
s = calc_elemfield(S,@sigmanodepc,q,varargin{:});
else
s = calc_elemfield(S,@sigmapc,q,varargin{:});
end


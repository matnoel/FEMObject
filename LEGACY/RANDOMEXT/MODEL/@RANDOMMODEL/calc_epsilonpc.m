function s=calc_epsilonpc(S,q,PC,varargin)

% function s=calc_epsilonpc(S,q,PC)
% S : MODEL
% q : deplacement
q=unfreevector(S,q);
if ischarin('node',varargin)
s = calc_elemfield(S,@epsilonnodepc,q,varargin{:});
else
s = calc_elemfield(S,@epsilonpc,q,varargin{:});
end


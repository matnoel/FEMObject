function f = calc_fintpc(S,q,PC,varargin)
% function f = calc_fintpc(S,q,PC,varargin)
% S : MODEL
% q : solution courante
% PC: POLYCHAOS

PC = getPC(PC);

q = unfreevector(S,q);
f = calc_pcvector(S,PC,@fintpc,q,varargin{:});
f = freevector(S,f);

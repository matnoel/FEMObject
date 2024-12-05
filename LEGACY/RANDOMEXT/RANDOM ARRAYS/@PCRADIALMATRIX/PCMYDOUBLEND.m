function u = PCMYDOUBLEND(u,k)

% function us = PCMYDOUBLEND(u,stodim)
% u : PCRADIALMATRIX
% us : PCMYDOUBLEND
% stodim : dimension stochastique du PCMYDOUBLEND (5 par defaut)

if nargin==1
    k=5;
end

a = MYDOUBLEND(u.V,k);
u = PCMYDOUBLEND(a,u.L,k);


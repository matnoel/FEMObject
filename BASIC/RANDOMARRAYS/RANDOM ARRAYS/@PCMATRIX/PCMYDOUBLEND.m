function u = PCMYDOUBLEND(u,k)

% function us = PCMYDOUBLEND(u,stodim)
% u : PCMATRIX
% us : PCMYDOUBLEND
% stodim : dimension stochastique du PCMYDOUBLEND (5 par defaut)


if nargin==1
    k=5;
end
if numel(u)==1
u = PCMYDOUBLEND(MYDOUBLEND(1),u,k);    
else
a = MYDOUBLEND(u.MULTIMATRIX,k);
u = PCMYDOUBLEND(a,getPC(u),k);
end

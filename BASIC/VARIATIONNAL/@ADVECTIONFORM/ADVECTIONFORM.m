function l = ADVECTIONFORM(c,tau,varargin)
% function a = ADVECTIONFORM(c,tau)
% a(v,u) = int tau v(c.grad(u))
% tau : champ par element ou champ par noeud ou double
%
if nargin==0
    l=struct();
    l.c = [];
    l.tau=1;
    l = class(l,'ADVECTIONFORM',setfree(BILINFORM(0,1),0));
else
    l=struct();
    l.c = c;
    if nargin==1
        l.tau = 1;
    else
        l.tau = tau;
    end
    l = class(l,'ADVECTIONFORM',setfree(BILINFORM(0,1),0));
end


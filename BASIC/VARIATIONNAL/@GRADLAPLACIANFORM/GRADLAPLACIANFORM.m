function l = GRADLAPLACIANFORM(tau,c1,c2,varargin)
% tau : champ par element ou champ par noeud ou double
% function a = GRADLAPLACIANFORM(tau,c1,c2)
% a(v,u) = int tau (c1.grad(v))(c2.delta(u))
%
% tau : champ par element ou champ par noeud ou double
%
if nargin==0
    l=struct();
    l.c1 = [];
    l.c2 = [];
    l.tau = 1;
    l = class(l,'GRADLAPLACIANFORM',setfree(BILINFORM(1,1),0));
else
    l=struct();
    l.tau = tau;
    if nargin==1
        l.c1 = [];
        l.c2 = [];
    elseif nargin==2
        l.c1 = c1;
        l.c2 = [];
    else
        l.c1 = c1;
        l.c2 = c2;
    end
    l = class(l,'GRADLAPLACIANFORM',setfree(BILINFORM(1,1),0));
end


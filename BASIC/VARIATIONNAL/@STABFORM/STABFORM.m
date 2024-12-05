function l = STABFORM(c,type,varargin)
% function a = STABFORM(c,type)
% stabilisation de la forme a(v,u) = int v(c.grad(u))
% tau : champ par element ou champ par noeud ou double
%
if nargin==0
    l=struct();
    l.c = [];
    l.tau=1;
    l.type=0;
    l.k = [];
    l = class(l,'STABFORM',setfree(BILINFORM(0,1),0));
else
    l=struct();
    l.c = c;
    l.tau = 1;
    l.type=type;
    if type==1
        l.k = varargin{1};
    else
        l.k=[];
    end

    l = class(l,'STABFORM',setfree(BILINFORM(0,1),0));
end


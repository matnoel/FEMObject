function l = MULTILINFORM(p,k,pk,varargin)
% function l = MULTILINFORM(p,k,pk)
% a(u1,u2,...) = int_domaine k . D^p(1)(u_1) D^p(2)(u_2) ...
%

if nargin==0
    l.p=[];
    l.k = [];
    l.pk=[];
    l = class(l,'MULTILINFORM',VARFORM());
else
    l=struct();
    l.p = p;
    if nargin<=1
        l.k = [];
    else
        l.k = k;
    end
    if nargin<=2
        l.pk=[];
    else
        l.pk=pk;
    end
    l = class(l,'MULTILINFORM',VARFORM(length(p),varargin{:}));
end


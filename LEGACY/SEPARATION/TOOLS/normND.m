function n = normND(u,varargin)
% function n = normND(u,varargin)

s=size(u);
if length(s)==2    
    if isa(u,'double') && nargin==1    
        n = norm(u,'fro');    
    else
        n = norm(u,varargin{:});
    end
else
    if nargin==1    
        n = u.*u;
    else
        n = u.*(varargin{1}*u);    
    end
    for k=length(s):-1:1
        n = sum(n,k);
    end
    n = sqrt(abs(squeeze(full(n))));
end

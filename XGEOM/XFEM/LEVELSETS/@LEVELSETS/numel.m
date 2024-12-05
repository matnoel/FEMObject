function s=numel(u,varargin)
if nargin==1
s=u.n;    
else
if isa(varargin{1},'char') && strcmp(varargin{1},':')
s=u.n;
else
s=length(varargin{1});    
end
end
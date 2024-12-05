function s=numel(u,varargin)
if nargin==1
s=u.M;    
else
if isa(varargin{1},'char') && strcmp(varargin{1},':')
s=u.M;
else
s=length(varargin{1});    
end
end
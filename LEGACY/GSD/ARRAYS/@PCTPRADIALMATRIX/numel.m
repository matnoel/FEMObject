function s=numel(u,varargin)
if nargin==1
s=prod(size(u));    
else
s=1; % pour quand on fait subsref avec {}
end
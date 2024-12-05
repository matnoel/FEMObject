function s=numel(u,varargin)
if nargin==1
    s=1;
else
    if isa(varargin{1},'char') && strcmp(varargin{1},':')
        s=1;
    else
        s=1;
    end
end
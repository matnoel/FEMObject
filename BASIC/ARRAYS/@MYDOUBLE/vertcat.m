function u = vertcat(varargin)
% function u = vertcat(varargin)

[rep,pos] = isclassin('MYDOUBLE',varargin);
u = varargin{pos(1)};
for i=1:nargin
    varargin{i} = double(varargin{i});
end
u.double = vertcat(varargin{:});

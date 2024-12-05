function a = VARFORM(n,varargin)

if nargin==0
    a.n = [];
    a.fact = [];
    a.selgroup = [];
    a.boundary = 0;
    a.boundaryobject = [];
    a.free = 1;
    a = class(a,'VARFORM');
else
    a.n = n;
    a.fact = 1;
    a.selgroup = getcharin('selgroup',varargin,'all');
    a.boundary = ischarin('boundary',varargin);
    if a.boundary
        a.boundaryobject = getcharin('boundary',varargin);
    else
        a.boundaryobject = [];
    end
    a.free = 1;
    a = class(a,'VARFORM');
end

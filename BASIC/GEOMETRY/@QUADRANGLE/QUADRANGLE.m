function u = QUADRANGLE(varargin)
% function D = QUADRANGLE(P1,P2,P3,P4)
% Pi : sommets du quadrangle

if nargin==0
    u = QUADRANGLE([0,0,0],[0,1,0],[1,1,0],[1,0,0]);
elseif nargin==1
    if isa(varargin{1},'QUADRANGLE')
        u = varargin{1};
    end
elseif nargin==4
    u.dim = 2;
    u.P{1} = POINT(varargin{1});
    u.P{2} = POINT(varargin{2});
    u.P{3} = POINT(varargin{3});
    u.P{4} = POINT(varargin{4});
    u.indim = getindim(u.P{1});
    u.diameter = max([distance(u.P{1},u.P{3}),...
        distance(u.P{1},u.P{4}),distance(u.P{1},u.P{2}),...
        distance(u.P{2},u.P{4})]);
    
    u = class(u,'QUADRANGLE',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end


function u = POLYGON(varargin)
% function D = POLYGON(P)
% P : argument POINT (plusieurs points)
% function D = POLYGON(P1,P2,...,Pn)
% P1,P2,...,Pn : argument POINT (plusieurs points)

if nargin==0
    u.dim = [];
    u.P = POINT();
    u.indim = getindim(POINT);
    
    u = class(u,'POLYGON',GEOMOBJECT(u.dim,u.indim));
    
elseif nargin==1
    if isa(varargin{1},'POLYGON')
        u = varargin{1};
    else
        u.dim = 2;
        u.P = POINT(varargin{1});
        u.indim = getindim(u.P);
        
        u = class(u,'POLYGON',GEOMOBJECT(u.dim,u.indim));
    end
    
else
    u.dim = 2;
    u.P = POINT(cell2mat(varargin(:)));
    u.indim = getindim(u.P);
    
    u = class(u,'POLYGON',GEOMOBJECT(u.dim,u.indim));
end

end


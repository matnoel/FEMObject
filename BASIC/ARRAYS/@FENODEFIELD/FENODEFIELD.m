function u = FENODEFIELD(varargin)
% constructeur de la classe FENODEFIELD
% u = FENODEFIELD(double)

if nargin==0
    u.size = [];
    u.value = [];
    
    u = class(u,'FENODEFIELD');
elseif isa(varargin{1},'FENODEFIELD')
    u = varargin{1};
else
    value = varargin{1};
    u.size = size(value);
    u.value = value;
    
    u = class(u,'FENODEFIELD');                    
end

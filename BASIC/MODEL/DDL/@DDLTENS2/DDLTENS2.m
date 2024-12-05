function v = DDLTENS2(varargin)
% function v = DDLTENS2(name,syscoord)
% name : 'EP' (deformation), 'SM' (contrainte), ...
%        'E' (deformation), 'G' (courbure), ...
%        'N' (effort), 'M' (moment), ...
% syscoord : objet de type systeme de coordonnee (CARTESIAN1D, CARTESIAN2D, ...)

switch nargin
    case 0
        v.name = [];
        v.ddl = DDL();
        v.nbddl = length(v.ddl);
        v.syscoord = SYSCOORD();
        
        v = class(v,'DDLTENS2');
        
    case 1
        if isa(varargin{1},'DDLTENS2')
            v=varargin{1};
        else
            v.name = varargin{1};
            v.ddl = DDL();
            v.nbddl = length(v.ddl);
            v.syscoord = SYSCOORD();
            
            v = class(v,'DDLTENS2');
        end
        
    case 2
        v.name = varargin{1};
        v.ddl = createddltens2(varargin{2},varargin{1});
        v.nbddl = length(v.ddl);
        v.syscoord = varargin{2};
        
        v = class(v,'DDLTENS2');
        
end

function v = DDLTENS4(varargin)
% function v = DDLTENS4(name,syscoord)
% name : 'C' (elasticite), 'S' (souplesse), ...
% syscoord : objet de type systeme de coordonnee (CARTESIAN1D, CARTESIAN2D, ...)

switch nargin
    case 0
        v.name = [];
        v.ddl = DDL();
        v.nbddl = length(v.ddl);
        v.syscoord = SYSCOORD();
        
        v = class(v,'DDLTENS4');
        
    case 1
        if isa(varargin{1},'DDLTENS4')
            v=varargin{1};
        else
            v.name = varargin{1};
            v.ddl = DDL();
            v.nbddl = length(v.ddl);
            v.syscoord = SYSCOORD();
            
            v = class(v,'DDLTENS4');
        end
        
    case 2
        v.name = varargin{1};
        v.ddl = createddltens4(varargin{2},varargin{1});
        v.nbddl = length(v.ddl);
        v.syscoord = varargin{2};
        
        v = class(v,'DDLTENS4');
        
end

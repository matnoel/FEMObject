function v = DDLVECT(varargin)
% function v = DDLVECT(name,syscoord,type)
% name : 'U' (deplacement), 'R' (rotation), ...
%        'F' (effort), 'M' (moment), ...
%        'EPS' (deformaton), 'GAM' (courbure), ...
%        'EFF' (effort), 'MOM' (moment), ...
%        'DU' (gradient du deplacement), 'DT' (gradient de temperature), ...
%        'Q' (flux), ...
% syscoord : objet de type systeme de coordonnee (CARTESIAN1D, CARTESIAN2D, ...)
% type : 'TRANS' ou 'ROTA'  (translation ou rotation)

switch nargin
    case 0
        v.name = '';
        v.type = '';
        v.ddl = '';
        v.nbddl = 0;
        v.syscoord = SYSCOORD();
        
        v = class(v,'DDLVECT');
        
    case 1
        if isa(varargin{1},'DDLVECT')
            v = varargin{1};
        else
            v.name = varargin{1};
            v.type = '';
            v.ddl = '';
            v.nbddl = 0;
            v.syscoord = SYSCOORD();
            
            v = class(v,'DDLVECT');
        end
        
    case {2,3}
        if nargin<=2
            type = 'TRANS';
        else
            type = varargin{3};
        end
        v.name = varargin{1};
        v.type = type;
        v.ddl = createddlvect(varargin{2},varargin{1},type);
        v.nbddl = length(v.ddl);
        v.syscoord = varargin{2};
        
        v = class(v,'DDLVECT');
        
    otherwise
        error('pas prevu')
        
end

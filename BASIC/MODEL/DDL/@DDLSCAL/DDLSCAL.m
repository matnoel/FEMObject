function v = DDLSCAL(varargin)
% function v = DDLSCAL(name)
% name : 'T' (temperature), 'U' (deplacement), 'QN' (deformation axiale ou flux normal), ...

switch nargin
    case 0
        v.ddl = {};
        v.nbddl = 0;
        v = class(v,'DDLSCAL');
        
    case 1
        if isa(varargin{1},'DDLSCAL')
            v = varargin{1};
        else
            v.ddl = {varargin{1}};
            v.nbddl = 1;
            v = class(v,'DDLSCAL');
        end
        
end

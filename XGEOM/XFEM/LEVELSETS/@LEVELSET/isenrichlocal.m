function e = isenrichlocal(c,varargin)

switch getnature(c)
    case 'material'
        e = getenrichtype(c)<=1;
    otherwise
        e=1;
end
function e = isenrichlocal(c,choix,varargin)

switch choix
    case {'support','cut'}
    e = getenrichtype(c,'support')~=1; 
    case {'tip','bicut'}
    e = getenrichtype(c,'tip',varargin{:})~=1;
    e = all(e);  
    otherwise
        error('pas defini')
end

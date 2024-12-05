function u = GMSHFILE(file)
% function u = GMSHFILE()

if nargin==0
    u.string = '';
    u.counter = 0;
    u.file = 'gmsh_file';
    u.iswritten = 0;
    u.ismesh = 0;
    
    u = class(u,'GMSHFILE');
elseif nargin==1 && isa(file,'GMSHFILE')
    u = file;
elseif nargin==1 && isa(file,'char')
    u = GMSHFILE();
    u = setfile(u,file);
end

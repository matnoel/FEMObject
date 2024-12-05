function u = setfile(u,file)
% function u = setfile(u,file)

if length(file)>=5 && strcmp(file(end-3:end),'.msh')
    u.ismesh = 1;
    u.file = file(1:end-4);
elseif length(file)>=6 && strcmp(file(end-4:end),'.mesh')
    u.ismesh = 1;
    u.file = file(1:end-5);
elseif length(file)>=5 && strcmp(file(end-3:end),'.geo')
    u.iswritten = 1;
    u.file = file(1:end-4);
else
    u.file = file;
end

function u = exportfile(u,options)
% function u = exportfile(u,options)

tempo = getfemobjectoptions('gmshpath');
commande = [tempo 'gmsh ' getfilemsh(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

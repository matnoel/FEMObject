function u = runfilemsh(u,options)
% function u = runfilemsh(u,options)

tempo = getfemobjectoptions('gmshpath');
commande = [tempo 'gmsh ' getfilemsh(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

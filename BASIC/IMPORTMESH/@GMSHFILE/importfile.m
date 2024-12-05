function u = importfile(u,options)
% function u = importfile(u,options)

tempo = getfemobjectoptions('gmshpath');
commande = [tempo 'gmsh ' getfilemesh(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

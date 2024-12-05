function u = runfile(u,options)
% function u = runfile(u,options)

tempo = getfemobjectoptions('gmshpath');
commande = [tempo 'gmsh ' getfilegeo(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

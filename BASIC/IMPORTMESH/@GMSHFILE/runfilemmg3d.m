function u = runfilemmg3d(u,options)
% function u = runfilemmg3d(u,options)

tempo = getfemobjectoptions('mmgpath');
commande = [tempo 'mmg3d_O3 ' getfilemsh(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

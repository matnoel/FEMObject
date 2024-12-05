function u = runfilemmg2d(u,options)
% function u = runfilemmg2d(u,options)

tempo = getfemobjectoptions('mmgpath');
commande = [tempo 'mmg2d_O3 ' getfilemsh(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

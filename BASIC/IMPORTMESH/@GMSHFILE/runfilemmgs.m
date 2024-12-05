function u = runfilemmgs(u,options)
% function u = runfilemmgs(u,options)

tempo = getfemobjectoptions('mmgpath');
commande = [tempo 'mmgs_O3 ' getfilemsh(u)];
if nargin==2
    commande = [commande ' ' options];
end
dos(commande);
u.ismesh = 1;

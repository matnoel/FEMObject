function initstate()
% function initstate()
% permet d'initialiser l'etat du generateur de nombres aleatoires

pause(eps)
rand('state',sum(100*clock));
randn('state',sum(100*clock));


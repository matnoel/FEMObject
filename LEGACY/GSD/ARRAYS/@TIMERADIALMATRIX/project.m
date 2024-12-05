function apc = project(apc,pc2)
% function apc = project(apc,pc2)
% projection de apc sur le POLYCHAOS pc2
% si apc est sur un chaos de dimension M1 et pc2 de dimension M2  
% si M1<M2  : repm localise pc2 dans pc1

apc.L = project(apc.L,pc2);
apc.POLYCHAOS = pc2;

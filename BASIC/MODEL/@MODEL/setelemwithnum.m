function [S,newgroup] = setelemwithnum(S,choix,info,numelem)
% function S = setlstypeelemwithnum(S,choix,type,numelem)
% pour chaque groupe d'element, on separe les elements dont le numero se 
% trouve dans la liste numelem et on modifie une des proprietes suivantes
% choix : 'lstype', 'lsenrich', 'lsenrichtype', 'lsnumber'
% info : valeur du parametre choix

[S,newgroup] = separateelemwithnum(S,numelem);
S = setgroupelemfield(S,choix,info,newgroup);



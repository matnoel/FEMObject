function [S,newgroup] = setparamelemwithnum(S,field,value,numelem)
% function S = setparamelemwithnum(S,choix,type,numelem)
% pour chaque groupe d'element, on separe les elements dont le numero se 
% trouve dans la liste numelem 
% et pour chaque nouveau groupe, on met la valeur value dans le parametre field


[S,newgroup] = separateelemwithnum(S,numelem);
S = setparamgroupelem(S,field, value,newgroup);

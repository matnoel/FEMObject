function [S,newgroup] = setfieldelemwithnum(S,field,value,numelem)
% function S = setfieldelemwithnum(S,choix,type,numelem)
% pour chaque groupe d'element, on separe les elements dont le numero se 
% trouve dans la liste numelem et on modifie une des proprietes suivantes
% et pour chaque nouveau groupe, on met la valeur value dans le champ field
%

[S,newgroup] = separateelemwithnum(S,numelem);
S = setfieldgroupelem(S,field, value,newgroup);

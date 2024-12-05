function dim2var=getdim2var(T)
%    function dim2var=getdim2var(T)
% Retourne le numero des variables contenues dans chaque
% dim de T, au premier niveau :

dim2var = cellfun(@(v2d) v2d(1),  T.var2dim ,'UniformOutput',1)';




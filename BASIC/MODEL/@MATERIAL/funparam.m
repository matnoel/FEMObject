function mat = funparam(mat,num,fun,varargin)
% function mat = funparam(mat,num,fun,varargin)
% appliquer la fonction fun aux parametres num du materiau
% num peut etre une liste de numero de parametres
% ou un char donnant le nom du parametre ou une cell de char
% varargin sont les arguments de fun

fun = fcnchk(fun);
if isa(num,'char')  || isa(num,'cell')
    [rep,num] = ischarin(num,mat.param(:,1));
    num = num(find(rep));
end
for i=1:length(num)
    mat.param{num(i),2} = fun(mat.param{num(i),2},varargin{:});
end


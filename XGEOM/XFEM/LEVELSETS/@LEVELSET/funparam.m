function ls = funparam(ls,num,fun,varargin)
% function ls = funparam(ls,num,fun,varargin)
% appliquer la fonction fun aux parametres num de la levelset
% num peut etre une liste de numero de parametres
% ou un char donnant le nom du parametre ou une cell de char
% varargin sont les arguments de fun

fun = fcnchk(fun);
if isa(num,'char')  || isa(num,'cell')
    [rep,num] = ischarin(num,ls.value(:,1));
    num = num(find(rep));
end
for i=1:length(num)
    ls.value{num(i),2} = fun(ls.value{num(i),2},varargin{:});
end


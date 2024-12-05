function [s,se] = integrate(S,funorder,fun,varargin)
% function [s,se] = integrate(S,funorder,fun,varargin)
% funorder : ordre d'integration ou function_handle telle que
% funorder(ELEMENT) rend l'ordre d'integration
%
% fun : fonction a integrer
% fun(xi,elem,xnode,varargin{:}) rend un MYDOUBLEND
% evaluation de la quantit� � int�grer en les points de Gauss xi

setemp = cell(1,getnbgroupelem(S));
node = S.node;

for p=1:getnbgroupelem(S)
    elem = getgroupelem(S,p);
    if isa(funorder,'function_handle')
        order = funorder(elem);
    else
        order = funorder;
    end
    xnode = node(elem);
    setemp{p} = integrate(elem,xnode,order,fun,varargin{:});
end

se = setemp{1};
for p=2:getnbgroupelem(S)
    se = concat(se,setemp{p},3);
end
s = sum(se,3);

if israndom(s)
    s = PCMATRIX(s);
else
    s = double(s);
end

function [s,se]=lsintegrate(S,funorder,funin,funout,varargin)
% function [s,se]=integrate(S,funorder,funin,funout,varargin)
% funorder : ordre d'integration ou function_handle telle que 
% funorder(ELEMENT) rend l'ordre d'integration
%
% funin et funout : fonction a integrer
% fun(xi,elem,xnode,varargin{:}) rend un MYDOUBLEND
% evaluation de la quantité à intégrer en les points de Gauss xi

setemp =cell(1,getnbgroupelem(S));
node=getnode(S);

for p=1:getnbgroupelem(S)
elem=getgroupelem(S,p);    
if isa(funorder,'function_handle')
    order = funorder(elem);
else
    order = funorder;
end
xnode = node(elem);
lsnum = getlsnumber(elem);
if ~isempty(lsnum)
 ls = getlevelset(S,lsnum);
 setemp{p} = lsintegrate(elem,xnode,order,ls,funin,funout,varargin{:});  
else
 setemp{p} = integrate(elem,xnode,order,funin,varargin{:});     
end

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




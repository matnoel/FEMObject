function isbicut = lsisbicut(elem,ls1,ls2,node)

if nargin==4
ls1xnode = getconnecvalue(elem,ls1,node);
ls2xnode = getconnecvalue(elem,ls2,node);
else
ls1xnode=ls1;
ls2xnode=ls2;
end

isbicut = lsiscut(elem,ls1xnode) & lsiscut(elem,ls2xnode);

if any(isbicut)
warning('peut ne pas etre bicut : à ameliorer')
end


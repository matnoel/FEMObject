function p = RANDPOLY(r)

if strcmp(func2str(r.X),'RVUNIFORM')
 p =POLYLEGENDRE();
elseif strcmp(func2str(r.X),'RVBETA')
if isa(r.funparam{1},'double') & isa(r.funparam{2},'double')
    p =POLYJACOBI(r.funparam{1:2});
else
    p=POLYLEGENDRE();
end

elseif strcmp(func2str(r.X),'RVGAMMA')
if isa(r.funparam{1},'double') 
    p =POLYLAGUERRE(r.funparam{1});
else
    p=POLYHERMITE();
end

else
    p=POLYHERMITE();
end

p=setnumber(p,getnumber(r));


function n = getn(a)

if length(a.forms)==0
    n=[];
else
n=getn(a.forms{1});
end



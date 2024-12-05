function s=size(u)
s=size(u.value);
if u.multidim>length(s)
 s = [s,ones(1,u.multidim-length(s))];
end    
s(u.multidim)=[];


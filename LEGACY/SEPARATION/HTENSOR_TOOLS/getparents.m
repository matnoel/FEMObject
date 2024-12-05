function p = getparents(u,t)
%function p = getparents(u,t)
p = t;
c = t;
while c~=1
    c = u.parent(c);
    p = [p c]; 
end
end

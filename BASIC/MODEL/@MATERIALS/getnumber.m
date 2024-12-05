function i = getnumber(u,k)
% function i = getnumber(u,k)

if nargin==1
    i = cell(0,u.n);
    for k=1:u.n
        i{k} = getnumber(u.MAT{k});
    end
else
    i = [];
    for j=1:length(k)
        i = [i,getnumber(u.MAT{k(j)})];
    end
end



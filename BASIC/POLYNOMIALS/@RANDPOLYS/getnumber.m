function i = getnumber(u,k)
% function i = getnumber(u,k)

if nargin==1
    i = cell(0,u.M);
    for k=1:u.M;
        i{k} = getnumber(u.h{k});
    end
else
    i = [];
    for j=1:length(k)
        i = [i,getnumber(u.h{k(j)})];
    end
end



function i = getnumber(u,k)
% function i = getnumber(u,k)

if nargin==1
    i = cell(0,u.n);
    for k=1:u.n;
        i{k} = getnumber(u.LS{k});
    end
else
    i = [];
    for j=1:length(k)
        i = [i,getnumber(u.LS{k(j)})];
    end
end



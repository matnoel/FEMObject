function s = size(u,dim)

if nargin==1
    s=zeros(1,u.dim);
    for i=1:u.dim
        if u.m(i)~=0
            s(i)=size(u.F{i}{1},1);
        else
            s(i)=0;
        end
    end
else
    s = size(u.F{dim}{1});
end

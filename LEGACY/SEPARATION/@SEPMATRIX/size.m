function s = size(u,dim)

if nargin==1
    s=zeros(1,u.dim);
    for i=1:u.dim
    s(i)=size(u.F{1,i},1);    
    end
else
    s = size(u.F{1,dim},1);
end





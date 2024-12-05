function s = issparse(u,dim)

if nargin==1
    s=zeros(1,u.dim);
    for i=1:u.dim
        s(i)=issparse(u.F{1,i});    
    end
else
    s = issparse(u.F{1,dim},1);
end


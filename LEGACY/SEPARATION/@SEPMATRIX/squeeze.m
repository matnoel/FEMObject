function u = squeeze(u)

s = size(u);
dim=find(s==1);
u.alpha = u.alpha.*prod(cell2mat(u.F(:,dim)),2)';
u.F(:,dim) = [];

if size(u.F,2)==0
    u = full(sum(u.alpha));  
elseif size(u.F,2)==1
    temp = u.alpha(1)*u.F{1};
    for k=2:size(u.F,1);
    temp = temp + u.alpha(k)*u.F{k};    
    end
    u = temp;
end

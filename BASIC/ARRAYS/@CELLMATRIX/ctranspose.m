function u=ctranspose(u)

if length(u.s)>=3
p=size(u.value,2);
s2 = u.s([2,1]);
s3 = u.s(3:end);
s = [s2,s3];
t = zeros(u.s(1:2));
t(:)=1:numel(t);
t=t';
t=t(:);
t = repmat(t(:),1,prod(s3));
j = repmat(1:prod(s3),prod(s2),1);
I = sub2ind([prod(s2),prod(s3)],t(:),j(:));

u.s=s;
u.value = conj(u.value(I,:));    
else
p=size(u.value,2);
s = u.s([2,1]);
t = zeros(u.s);
t(:)=1:numel(t);
t=t';

I = ind2sub(s,t(:));
u.s=s;

u.value = conj(u.value(I,:));
end
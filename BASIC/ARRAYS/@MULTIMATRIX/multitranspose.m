function u=multitranspose(u)
%warning('attention c''est redefini')
if isa(u.value,'cell')
   u.value = u.value'; 
   u.sm = size(u.value);
else
sm = u.sm([2,1]);
t = zeros(u.sm);
t(:)=1:numel(t);
t=t';

I = ind2sub(sm,t(:));
u.sm=sm;

u.value = u.value(:,I);

end
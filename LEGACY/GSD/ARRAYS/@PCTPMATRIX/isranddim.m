function a = isranddim(x,i)


if nargin==1
a = x.isranddim;
else
a = x.isranddim(i);    
end

return

if nargin==1
i = 1:getnbdim(x);
end
a = zeros(1,length(i));
for k=1:length(i)
    a(k) = numel(x.phi{i(k)})>1;
end

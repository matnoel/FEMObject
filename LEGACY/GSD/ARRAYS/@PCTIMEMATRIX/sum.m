function v = sum(u,k)

if nargin==2
    if k~=1
        error('pas programme')
    end
    v=u;
    v.value = sum(u.value,k);
    v.s = size(v.value);
    v.s(2)=1;
else
v = sum(u,1);
end
    
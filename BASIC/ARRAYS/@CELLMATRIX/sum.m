function u=sum(u,k)

if nargin==1 & any(u.s==1)
    k=0;
else
    k=1;
end

s=u.s;
sm=u.sm;
switch k
    case 0
    u.value = sum(u.value,1);
    u.s = [1,1]; 
    case 1
     u.value = sum(reshape(u.value,s(1),s(2)*prod(sm),1)); 
     u.s = [1,s(2)];
     u.value = reshape(u.value,s(2),prod(sm));
    case 2
     u.value = sum(reshape(u.value.',s(1)*prod(sm),s(2),2)); 
     u.s = [s(1),1];
     u.value = reshape(u.value,s(1),prod(sm)).';     
    otherwise
    error('utiliser multisum')
end

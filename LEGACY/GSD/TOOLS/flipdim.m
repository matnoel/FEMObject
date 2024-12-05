function u = flipdim(u,k)

n = size(u,k);
f = n:-1:1;

switch k
case 1
    u = u(f,:);  
case 2
    u = u(:,f);  
otherwise
    error('pas prevu')
end


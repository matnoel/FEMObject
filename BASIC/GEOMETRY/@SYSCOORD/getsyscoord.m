function v = getsyscoord(v,p)
% function v = getsyscoord(v,p)

if nargin==2 && size(v.MYDOUBLEND,3)>1    
    v.MYDOUBLEND = v.MYDOUBLEND(:,:,p);
end

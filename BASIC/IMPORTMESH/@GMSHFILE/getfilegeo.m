function f = getfilegeo(u,iter)
% function u = getfilegeo(u,iter)

if nargin==1 || isempty(iter)
    f = [u.file '.geo'];
else
    f = [u.file '_' num2str(iter) '.geo'];
end

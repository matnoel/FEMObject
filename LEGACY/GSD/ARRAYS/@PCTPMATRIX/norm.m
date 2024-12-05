function s = norm(x)
% function s = norm(x)

s = norm(x.phi0);
for i=1:getnbdim(x)
    s = s*norm(x.phi{i});
end

function e = getnature(ls)
% function e = getnature(ls)

if isempty(ls.material)
    e = 'domain';
else
    e = 'material';
end
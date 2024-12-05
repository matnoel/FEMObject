function u=double(u)
if iseval(u)
    u=double(u.value);
else
    error('u must be a LEVELSET')
end
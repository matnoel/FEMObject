function disp(a)
% function disp(a)

s = ['a(u,v) = int( ' getstring('k',a.pk) '.' getstring('u',a.p) '.' getstring('v',a.q) ' )'] ;
disp(s);
disp(' ')


function s = getstring(l,p)
% function s = getstring(l,p)

if isempty(p) || p==0
    s = l;
elseif p==1
    s = ['D(' l ')'];
else
    s = ['D^' num2str(p) '(' l ')'];
end

return

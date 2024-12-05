function display(u)
% function display(u)

if length(u.funs)==0
    s = [0,0];
else
    s = size(u);
end
siz = num2str(s(1));
for k=2:length(s)
    siz = [siz '-by-' num2str(s(k))];
end

disp(' ')
disp([inputname(1) ' = ( ' siz ' ' class(u) ' )'])
disp(' ')
disp(u)

function display(u)
% function display(u)

s = size(u);
siz = num2str(s(1));
for k=2:length(s)
    siz = [siz '-by-' num2str(s(k))];
end

if isa(u.value,'MULTIMATRIX')
    addinfo = ' (MULTIMATRIX) ';
else
    addinfo = ' ';
end

disp(' ')
disp([inputname(1) ' = ( ' siz ' ' class(u) ' )' addinfo])
disp(' ')
disp(u)


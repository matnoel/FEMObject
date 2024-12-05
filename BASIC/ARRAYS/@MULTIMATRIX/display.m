function display(u,varargin)
% function display(u,varargin)

siz = num2str(u.s(1));
for k=2:length(u.s)
    siz = [siz '-by-' num2str(u.s(k))];
end
sizm = num2str(u.sm(1));
for k=2:length(u.sm)
    sizm = [sizm '-by-' num2str(u.sm(k))];
end

disp(' ')
disp([inputname(1) ' = ( ' siz ' ' class(u) ' ' sizm ' ) '])
disp(' ')

if ~ischarin('novalue',varargin)
    disp(u)
end
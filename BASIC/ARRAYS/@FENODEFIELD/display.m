function display(u)
% function display(u)

disp(' ')
disp([inputname(1) ' = (FENODEFIELD)'])
if isa(u.value,'double')
    disp(' ')
    disp(u.value)
else
    value = u.value;
    display(value)
end

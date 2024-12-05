function elem = setsyscoord(elem,syscoord)
% function elem = setsyscoord(elem,syscoord)
%
if isa(syscoord,'SYSCOORD')
    elem.syscoord = syscoord;
else
    error(' ')
end






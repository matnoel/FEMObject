function elem = setsyscoordlocal(elem,syscoord)
% function elem = setsyscoordlocal(elem,syscoord)
% 
if isa(syscoord,'SYSCOORD')
    elem.syscoordlocal = syscoord;
else
    error(' ')
end






function s = entity(name,number,values)
% function s = entity(name,number,values)

s = [name '(' num2str(number) ') = '];

s = [s, valuesintobraces(values)];

s = [s, ';\n'];

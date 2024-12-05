function u = setfieldattribute(u,name,value,number)
% function u = setfieldattribute(u,name,value,number)

if nargin<4 || isempty(number)
    number = 1;
end

u = stepcounter(u);
s = ['Field[' num2str(number) '].' name ' = ' num2str(value) ' ;\n'];
u = addstring(u,s);

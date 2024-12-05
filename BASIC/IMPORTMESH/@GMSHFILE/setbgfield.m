function u = setbgfield(u,number)
% function u = setbgfield(u,number)

if nargin<4 || isempty(number)
    number = 1;
end

u = stepcounter(u);
s = ['Background Field = ' num2str(number) ' ;\n'];
u = addstring(u,s);

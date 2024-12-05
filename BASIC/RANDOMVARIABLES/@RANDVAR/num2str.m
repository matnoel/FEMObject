function s=num2str(u)

s=['RANDVAR #' num2str(getnumber(u)) ' -> ' class(u) '('];

for i=1:size(u.param,1)
    if i>1
        s=[s,','];
    end
       s=[s,num2str(u.param{i,2})]; 
end


s=[s ')'];
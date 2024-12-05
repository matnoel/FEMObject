function s = valuesintobraces(values)
% function s = valuesintobraces(values)

s = [' {' num2str(values(1))];
for k=2:length(values)
    s = [ s, ' , '  num2str(values(k)) ];
end
s = [s, '} '];

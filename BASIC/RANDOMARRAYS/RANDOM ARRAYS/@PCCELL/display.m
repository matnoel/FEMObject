function display(upc)
% function display(upc)

disp(' ')
disp([inputname(1) ' = (' class(upc) ')'])
disp(' ')
for i=0:getP(upc)
    fprintf([ inputname(1) '{%d} = \n\n' ],i)
    disp(upc.value{i+1});
    fprintf('\n')
end

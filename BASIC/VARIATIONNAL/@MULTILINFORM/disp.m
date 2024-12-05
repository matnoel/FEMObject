function disp(a)
% function disp(a)

s = inputname(1);
switch getn(a)
    case 1
        sarg = {'v'};
    case 2
        sarg = {'v','u'};
    case 3
        sarg = {'u','v','w'};
    otherwise
        sarg = {};
        for i=1:getn(a)
            sarg = [sarg , {['u' num2str(i)]}];
        end
end

s = [inputname(1) '(' sarg{1}];
for i=2:getn(a)
    s = [s ',' sarg{i}];
end

if ~isempty(a.k)
    s = [s ') = int( k.'];
else
    s = [s ') = int( '];
end


s = [s getstring(sarg{1},a.p(1))];
for i=2:getn(a)
    s = [s '.' getstring(sarg{i},a.p(i))];
end

s = [s ' )'];
disp(s)
disp(' ')


function s = getstring(l,p)
% function s = getstring(l,p)

if isempty(p) || p==0
    s = l;
elseif p==1
    s = ['grad(' l ')'];
else
    s = ['D^' num2str(p) '(' l ')'];
end

return

function v = subsref(u,s)
% function v = subsref(u,s)

comm = 'u';

for l=1:length(s)
    
    switch s(l).type
        case '.'
            comm = [comm '.' s(l).subs];
        case '{}'
            comm = [comm '{s(' num2str(l) ').subs{1}'];
            for i=2:length(s(l).subs)
                comm = [comm ',s(' num2str(l) ').subs{' num2str(i) '}'];
            end
            comm = [comm '}'];
        case '()'
            comm = [comm '(s(' num2str(l) ').subs{1}'];
            for i=2:length(s(l).subs)
                comm = [comm ',s(' num2str(l) ').subs{' num2str(i) '}'];
            end
            comm = [comm ')'];
    end
end

v = eval(comm);

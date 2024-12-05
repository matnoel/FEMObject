function disp(u)
% function disp(u)
% Peut etre amelioree :
%   afficher la dimension, la nature des variables etc

if u.dim~=0
    for d=1:u.dim-1
        if isa(u.F{d}.model,'double')
            aff = ['[COLOCATION]    '];
        else
            aff = ['[' class(u.F{d}.model) ']    '];
        end
        fprintf(aff);
    end
    d = u.dim;
    disp(['[' class(u.F{d}.model) ']']);
    fprintf('\n')
else
    disp('SEPMODEL vide');
end


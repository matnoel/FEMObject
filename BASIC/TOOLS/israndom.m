function x = israndom(u)
% function x = israndom(u)
% x = 1 si u est aleatoire
%     0 sinon

if isa(u,'cell')
    x = false;
    for i=1:length(u)
        x = x | israndom(u{i});
    end
else
    x = 0;
end

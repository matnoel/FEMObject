function r=cellstrcmp(s1,s2)
% function r=cellstrcmp(s1,s2)
% r=1 si tous les couples de cellules (s1{i}, s2{i}) contiennent les memes chaines de caracteres
% r=0 sinon (notamment s'il n'y a pas le meme nombre de cellules dans s1 et s2)
r=1;
if length(s1)~=length(s2)
    r=0;
else
    for i=1:length(s1)
        if ~strcmp(s1{i},s2{i})
            r=0;
            break
        end
    end
end


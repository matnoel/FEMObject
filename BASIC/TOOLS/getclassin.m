function [v,pos,num] = getclassin(s,var,default)
% function [v,pos] = getclassin(classname,var,default)
% var : tableau de cellules
% v : objet de classe classname (ou heritant de classname) contenu dans var
%     si plusieurs objets, v est une cellule
% pos : position de l'objet (ou des objets) dans le tableau var

pos=[];
rep=0;
for i=1:length(var)
    if isa(var{i},s)
        rep=1;
        pos=[pos,i];
    end
end
if length(pos)==1
    v=var{pos};
elseif length(pos)>1
    v=var(pos);
else
    v=[];
end

num = length(pos);

if isempty(pos) && nargin>2
    v=default;
end

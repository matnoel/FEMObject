function s=copyfields(s,a)
% function s=copyfields(s,a)
% s : structure 
% a : structure ou objet (dans ce cas il est converti en structure)
% copie les champs de a dans la structure s


a=struct(a);
f=[fieldnames(a)';struct2cell(a)'];
s=setfield(s,f{:});



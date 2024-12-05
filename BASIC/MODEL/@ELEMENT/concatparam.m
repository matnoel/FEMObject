function elem1 = concatparam(elem1,elem2)
% concatene les parametres des elements

param1 = getparam(elem1);
param2 = getparam(elem2);
f1=fieldnames(param1);
f2=fieldnames(param2);
param1 = struct2cell(param1);
param2 = struct2cell(param2);
param = [f1(:),param1(:)]';

if length(f1)~=length(f2)
    error('il doit y aoir le meme nombre de parametres')
end

for i=1:length(f1)
    try
        param{2,i}=concat(param1{i},param2{i});
%     catch
%         warning('un parametre ne peut etre concatener');
    end
end

param=param;

param = struct(param{:});
elem1=setparam(elem1,param);

if ~isempty(elem1.numddl) & ~isempty(elem2.numddl)
    elem1.numddl=[elem1.numddl;elem2.numddl];
end


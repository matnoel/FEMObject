function S = setgroupelemlstype(S,lstype,num)
% function S = setlstypeelemwithnum(S,lstype,numgroup)
% associe le 'lstype' lstype aux groupes d'elements numgroup

if nargin==2
num = 1:length(S.groupelem);
end
if ~isempty(num) && max(num)>length(S.groupelem)
    error('les groupes d''elements ne''existent pas')
end
for i=1:length(num)
 S.groupelem{num(i)} =  setlstype(S.groupelem{num(i)},lstype);
end



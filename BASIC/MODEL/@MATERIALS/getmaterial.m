function mat = getmaterial(mats,k)
% function mat = getmaterial(mats,k)

if nargin==1
    mat = mats.MAT;
else
    [ok,rep] = ismember(k,mats);
    if length(k)~=length(rep) || any(rep==0)
        error('Wrong number')
    end
    mat = mats.MAT(rep);
end

if length(mat)==1
    mat = mat{1};
elseif isempty(mat)
    mat=[];
end

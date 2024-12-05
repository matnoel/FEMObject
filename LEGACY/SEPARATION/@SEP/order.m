function o=order(A)
% function o=order(A)
% o est l'ordre du tenseur A
% Si o=-1, les ordres ne st pas les memes
o=zeros(1,A.dim);
for d=1:A.dim
    if isa(A.F{1,d},'double')
        o(d)=length(size(A.F{1,d}));
    elseif isa(A.F{1,d},'MULTIMATRIX')
        o(d)=3; % A priori...
    elseif isa(A.F{1,d},'SPARSETENSOR')
        o(d)=A.F{1,d}.ordre;
    end
end
o=unique(o);
if length(o)~=1
    o=-1;
end
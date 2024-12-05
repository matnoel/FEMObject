function s = size(H)
s=[];
for i=1:H.dim
    s=[s size(H.F{1,i})];
end


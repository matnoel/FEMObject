function s = issparse(H)
s=[];
for i=1:H.dim
    s=[s issparse(H.F{1,i})];
end


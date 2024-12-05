function H=hseprand(T,s)

S = SEPMATRIX(length(s),1);
for d=1:length(s)
    S(d)=rand(s(d),1);
end

H = HSEPMATRIX(S,T);
% H = normalize(H);

function r=gettypes(B)

r=zeros(1,length(B.BC));

for i=1:length(B.BC)
    r(i)= B.BC{i}.type;
end
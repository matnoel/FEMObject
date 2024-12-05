function a=mean(pcr)

z=mean(vertcat(pcr.L{:}));
z=pcr.D*z;

a=pcr.V{1}*z(1);
for i=2:pcr.m
a=a+ pcr.V{i}*z(i);   
end



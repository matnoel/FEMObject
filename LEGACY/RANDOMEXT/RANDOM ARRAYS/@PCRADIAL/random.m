function a = random(pcr)
% function a = random(pcr)

z = pcr.D*random(vertcat(pcr.L{:}),1);

a = pcr.V{1}*z(1);
for i=2:pcr.m
    a = a+pcr.V{i}*z(i);
end


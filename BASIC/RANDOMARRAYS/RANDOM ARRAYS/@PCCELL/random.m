function a = random(apc,varargin)
% function a = random(apc,varargin)

rep = ischarin('init',varargin);
if rep
    initstate
end
D = random(apc.POLYCHAOS,n);
a = D(1)*apc.value{1};
for i=1:getP(PC)
    a = a+D(i+1)*apc.value{i+1};
end

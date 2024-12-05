function [a,xi]=random(A,varargin)
% function [a,xi]=random(A,varargin)

n = [varargin{:}];
Lb = A.L;
[Lb{1},xi] = random(Lb{1},prod(n));

for i=2:length(Lb)
    Lb{i} = randomeval(Lb{i},xi);
end

a = A.V{1}*Lb{1};
for k=2:A.m
    a = a + A.V{k}*Lb{k};
end

function [ok,rep] = israndomparam(param,k)
% function [ok,rep] = israndomparam(param,k)

if nargin==1
    k=1:size(param,1);
end

rep=zeros(1,length(k));
for i=1:length(k)
    rep(i) = israndom(param{k(i),2});
end

ok = any(rep);
rep = find(rep);

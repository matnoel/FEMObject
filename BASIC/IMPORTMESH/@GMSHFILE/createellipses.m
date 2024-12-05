function u = createellipses(u,numbercenter,seg,maj,numberlines)
% function u = createellipses(u,numbercenter,seg,maj,numberlines)

for k=1:size(seg,1)
    u = createellipse(u,numbercenter,seg(k,:),maj(k),numberlines(k));
end

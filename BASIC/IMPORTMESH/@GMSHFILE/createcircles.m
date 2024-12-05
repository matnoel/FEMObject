function u = createcircles(u,numbercenter,seg,numberlines)
% function u = createcircles(u,numbercenter,seg,numberlines)

for k=1:size(seg,1)
    u = createcircle(u,numbercenter,seg(k,:),numberlines(k));
end

function u = createlines(u,seg,numberlines)
% function u = createlines(u,seg,numberlines)

for k=1:size(seg,1)
    u = createline(u,seg(k,:),numberlines(k));
end

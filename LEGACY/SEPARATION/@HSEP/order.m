function o=order(H)
o=cellfun(@(h) order(h) ,H.F(1,:));
o=unique(o);
if length(o)~=1
    o=-1;
end
function u = sptensordiag(v,dim)
% function u = tensordiag(v)

if ~isa(v,'double') || ~numel(find(size(v)>1))>1
    error('enter a vector')
end
m = length(v);
m = repmat(m,1,dim);
u= sptensor(m);
s.type = '()';       
for i=1:length(v)                         
    s.subs = repmat({i},1,dim);
    u = subsasgn(u,s,v(i));
end





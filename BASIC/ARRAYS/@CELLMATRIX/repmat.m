function w=repmat(w,a,b)

if nargin==2
if length(a)~=2
    error('mauvais argument dans repmat')
end
b=a(2);
a=a(1);
end

m=size(w.value,2);
if a>1
w.value = reshape(w.value,w.s(1),w.s(2)*m);
w.value = repmat(w.value,a,1);
w.s(1) = size(w.value,1);
w.value = reshape(w.value,w.s(1)*w.s(2),m);
end
if b>1
w.value = w.value.';
w.value = reshape(w.value,m*w.s(1),w.s(2));
w.value = repmat(w.value,1,b);
w.s(2) = size(w.value,2);
w.value = reshape(w.value,m,w.s(1)*w.s(2)).';
end
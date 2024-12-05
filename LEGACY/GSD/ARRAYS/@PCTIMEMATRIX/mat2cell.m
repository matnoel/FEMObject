function w = mat2cell(u)

w=u;    
if isa(u.value,'cell')
    warning('PCTIMEMATRIX deja stocke sous forme de cell')
else
n= length(u.TIMEMODEL);
w.value = cell(1,n);
s  = u.s ; 
if s(2)>1
   error('a verifier') 
end
for k=1:n
w.value{k} = reshape(u.value(:,k),s);        
end
end
function a=setfree(a,i)
if ~(i==0 || i==1)
    error('mettre 1 ou 0')
end
a.free = i;
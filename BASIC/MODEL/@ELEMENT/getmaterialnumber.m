function n=getmaterialnumber(elem)

if isempty(elem.material)
n=0;    
else
n=getnumber(elem.material);
end


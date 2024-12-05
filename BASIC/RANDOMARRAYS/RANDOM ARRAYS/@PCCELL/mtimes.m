function w=mtimes(u,v)

if isa(u,'PCCELL') & ~isa(v,'PCCELL')
   w=u;
   
   if isa(v,'PCARRAY')
   s=size(w);
   v=double(v);
   masse=getmasse(u);
    for k=1:length(masse)
    massev{k} = masse{k}*v';
    end
    massev = getmatrix(massev);
    w.value = getmatrix(w.value);
    w.value = w.value*massev';
    w.value = setmatrix(w.value,s(1),s(2));
   else
        
        for k=1:length(w.value)
         w.value{k} = w.value{k} * v ;   
        end
   
   end
    
elseif isa(v,'PCCELL') & ~isa(u,'PCCELL')
    w=v;
    if isa(u,'PCARRAY')
   s=size(w);
   u=double(u);
   masse=getmasse(v);
    for k=1:length(masse)
    massev{k} = masse{k}*u';
    end
    massev = getmatrix(massev);
    w.value = getmatrix(w.value);
    w.value = w.value*massev';
    w.value = setmatrix(w.value,s(1),s(2));
    else
      
        for k=1:length(w.value)
         w.value{k} = u*w.value{k} ;   
        end
   
    end
    
else
       error('mtimes not defined') 
end


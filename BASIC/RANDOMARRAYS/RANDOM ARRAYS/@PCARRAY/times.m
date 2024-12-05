function w=times(u,v)

if isa(u,'PCARRAY') & isa(v,'double') 
    w=PCRADIAL(v,u);
elseif isa(v,'PCARRAY') & isa(u,'double') 
    w=PCRADIAL(u,v);
else
    error('times not defined')
end

function xnode=calc_xnode(elem,xnode)



if isa(xnode,'NODE')
    xnode = getcoord(xnode,getconnec(elem)');
end
if isa(xnode,'double')
    xnode= MYDOUBLEND(xnode);
end


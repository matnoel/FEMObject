function x = calc_midpointelem(elem,node)
[temp,connec] = ismember(getconnec(elem),getnumber(node)) ;
x=mean(getcoord(node,elem),1);
x=POINT(x);
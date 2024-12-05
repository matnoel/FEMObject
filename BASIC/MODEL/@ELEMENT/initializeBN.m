function elem = initializeBN(elem,node,varargin)

n = getcharin('intorder',varargin,'rigi');
elem = setparam(elem,'intorder',n);
xnode = node(elem);
gauss=calc_gauss(elem,n);
elem = setparam(elem,'gauss',gauss);
[B,detJ]=calc_B(elem,xnode,gauss.coord);
[N,detJ,x]=calc_N(elem,xnode,gauss.coord);
elem = setparam(elem,'B',B);
elem = setparam(elem,'N',N);
elem = setparam(elem,'detJ',detJ);
elem = setparam(elem,'x',x);

elem = setparam(elem,'initializeBN',true);


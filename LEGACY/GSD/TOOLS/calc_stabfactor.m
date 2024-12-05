function tau = calc_stabfactor(S,k,c)

node = getnode(S);
tau = zeros(getnbelem(S),1);
c = full(c);
for i=1:getnbgroupelem(S)
    elem=getgroupelem(S,i);
    xnode = node(elem);
    gauss=calc_gauss(elem,0);
    B = calc_B(elem,xnode,gauss.coord);
    if numel(k)==getnbnode(S) || numel(c)==getnbnode(S)*getdim(elem)
        N = calc_N(elem,xnode,gauss.coord);
    end

    if numel(k)==1
        ke=k;
    elseif numel(k)==getnbnode(S)
        ke = localize(elem,k);
        ke = N*full(ke);
    end

    if numel(c)==getdim(elem)
        ce=c;
        bxi = norm(c);
    elseif numel(c)==getnbnode(S)*getdim(elem)
        ce = localize(elem,c);
        ce = N*ce;
        bxi = norm(ce);
    else
        keyboard
    end

    he = 2/sum(abs(ce'*B/bxi));
    Pe = (bxi*he/2*ke^-1); 
    taue = he/2./bxi.*(1./tanh(Pe)-1./Pe);

    tau(getnumber(elem))=taue;

end

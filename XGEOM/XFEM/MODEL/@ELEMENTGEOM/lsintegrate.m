function I = lsintegrate(elem,xnode,order,ls,funin,funout,varargin)

switch getlstype(elem)
    case 'in'
      I = integrate(elem,xnode,order,funin,varargin{:});         
    case 'out'
      I = integrate(elem,xnode,order,funout,varargin{:});           
    case 'cut'
connec = getconnec(elem);
lsval = getvalue(ls);
PC=[];
for e=1:getnbelem(elem)
    eleme = getelem(elem,e);
    connece = connec(e,:);
    xnodee = xnode(:,:,e);
    lse = lsval(connece);
    
    if all(lse<=0)
        gauss = calc_gauss(eleme,order);
        Ie = integrate_with_gauss(eleme,xnodee,gauss,funin,varargin{:});
    elseif all(lse>=0)
        gauss = calc_gauss(eleme,order);
        Ie = integrate_with_gauss(eleme,xnodee,gauss,funout,varargin{:});
    else
        [gaussin,gaussout] = calc_lssubgauss(eleme,lse,order);
        Ie = integrate_with_gauss(eleme,xnodee,gaussin,funin,varargin{:}) + ...
             integrate_with_gauss(eleme,xnodee,gaussout,funout,varargin{:});
    end
    if isa(Ie,'PCMYDOUBLEND')
        Ie=expand(Ie);
        Ie = permutestodim(Ie,getstodim(Ie),5);
    end
    if e==1 && ~israndom(Ie)
        I = zerosND([size(Ie),getnbelem(elem)]);
    elseif e==1 && isa(Ie,'PCMYDOUBLEND')      
        I = zerosND([size(Ie),getnbelem(elem),1,length(getPC(Ie))]);
    end
    
    if ~israndom(Ie)
        I(:,:,e) = Ie;    
    else
        Ie = getV(expand(Ie));
        PC = getPC(Ie);
        for k=1:length(PC)
        I(:,:,e,1,k) = Ie(:,:,1,1,k);
        end
    end
    
end
    if ~isempty(PC)
        I = PCMYDOUBLEND(I,PC,5);
    end
end

if isa(I,'PCMYDOUBLEND')
    I = permutestodim(I,getstodim(I),5);
end


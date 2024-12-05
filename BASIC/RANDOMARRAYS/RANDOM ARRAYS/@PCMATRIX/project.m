function apc2 = project(apc,pc2)
% function apc2 = project(apc,pc2)
% projection de apc sur le POLYCHAOS pc2
% si apc est sur un chaos de dimension M1 et pc2 de dimension M2
% si M1<M2  : repm localise pc2 dans pc1

pc1 = apc.POLYCHAOS;
pc2 = getPC(pc2);
H1=RANDPOLYS(pc1);
H2=RANDPOLYS(pc2);
if getM(pc1)<=getM(pc2)
    [ok,repm] = ismember(H1,H2); % repm localise pc1 dans pc2
    repm = repm(find(ok));
    
else
    [ok,repm] = ismember(H2,H1);% repm localise pc2 dans pc1
    repm = repm(find(ok));
end


if (getM(pc1)<=getM(pc2) && ~polycmp(H1,H2(repm))) || (getM(pc1)>getM(pc2) && ~polycmp(H1(repm),H2))
    
    if getM(pc1)<=getM(pc2)
        % fun = inline('MULTIMATRIX(randomeval(fpc,transfer(x1,x2(repm),a(:,repm))),size(fpc))','a','fpc','x1','x2','repm');
        fun = @(a,fpc,x1,x2,repm) MULTIMATRIX(randomeval(fpc,transfer(x1,x2(repm),a(:,repm))),size(fpc));
    else
        % fun = inline('MULTIMATRIX(randomeval(fpc,transfer(x1(repm),x2,a)),size(fpc))','a','fpc','x1','x2','repm');
        fun = @(a,fpc,x1,x2,repm) MULTIMATRIX(randomeval(fpc,transfer(x1(repm),x2,a)),size(fpc));
    end
    Hint = intersect(H2,H1);
    ng = round((1+max(getorder(pc1))+max(getorder(pc2)))/2)+1;
    apc2 = decompmatrix(pc2,ng,Hint,fun,apc,RANDVARS(pc2),RANDVARS(pc1),repm);
    
else
    apc1 = double(apc);
    apc2 = sparse(size(apc1,1),length(pc2));
    
    ind1=getindices(pc1);ind1=ind1(:,1:end-1);
    ind2=getindices(pc2);ind2=ind2(:,1:end-1);
    
    if getM(pc1)==getM(pc2)
        
        [int,i1,i2]=intersect(ind1,ind2,'rows');
        apc2(:,i2) = apc1(:,i1);
        
        apc2 = PCMATRIX(apc2,size(apc),pc2);
        
    elseif getM(pc1)<getM(pc2)
        repnm=setdiff(1:getM(pc2),repm);
        %[temp,indices]=ismember(ind2(:,repm),ind1,'rows');
        %rep = find(temp);
        pc2moins1 = restrictdim(pc2,getnumber(pc2,repnm));
        pc2moins1 = setindices(pc2moins1,ind2(:,repnm));
        Hmprod = mean(pc2moins1);
        %moy2 = double(one(pc2moins1));
        %H = RANDPOLYS(pc2);
        %Hm=sparse(getM(H),size(ind2,1));
        %for k=repnm
        %    Hm(k,:) = mean(H{k},ind2(:,k))';
        %end
        %Hmprod = prod(Hm(repnm,:),1);
        %moy = double(one(pc2));
        %apc2(:,rep)=apc1(:,indices(rep));
        %apc2 = apc2 .* repmat(Hmprod,[size(apc2,1),1]);
        
        Hmprod = Hmprod(:);
        repH = find(Hmprod);
        [temp,indices]=ismember(ind2(repH,repm),ind1,'rows');
        rep = find(temp);
        
        apc2(:,repH(rep))=apc1(:,indices(rep)).*repmat(Hmprod(repH(rep))',size(apc1,1),1);
        
        apc2 = PCMATRIX(apc2,size(apc),pc2);
        
    else
        H = RANDPOLYS(pc1);
        repnm=setdiff(1:getM(pc2),repm);
        Hm=sparse(getM(H),size(ind1,1));
        for k=repnm
            Hm(k,:) = mean(H{k},ind1(:,k))';
        end
        Hm = prod(Hm(repnm,:),1);
        for j=1:size(ind2,1)
            [temp,indices]=ismember(ind1(:,repm),ind2(j,:),'rows');
            rep = find(temp);
            apc2(:,j) = apc1(:,rep)*Hm(rep)';
        end
        apc2 = PCMATRIX(apc2,size(apc),pc2);
        
    end
    
end

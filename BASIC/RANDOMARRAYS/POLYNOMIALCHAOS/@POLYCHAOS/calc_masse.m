function PC=calc_masse(PC,PC2)
% function PC=calc_masse(PC,PC2)
% PC, PC2 : POLYCHAOS
% masse{alpha} matrice de taille (P2+1)-by-(P2+1)   (alpha=1...P+1)
% masse{alpha}(beta,gamma)=E(H_alpha H_beta H_gamma)


if nargin==1
    PC2=getPC(PC);
elseif ~isa(PC2,'POLYCHAOS')
error('arguments POLYCHAOS')  
end

try
  % load('MULTIMASSE') 
end
ismasse = 0;isin1 = 0; isin2=0;
if exist('masse')
ismasse = 1;  
try
isin1 = isin(PC,PCsave);
isin2 = isin(PC2,PC2save);
end
end

if ismasse && isin1 && isin2
%fprintf('-> Chargement de la matrice masse stochastique\n')

[rep,ia]=isin(PC,PCsave);
[rep2,ia2]=isin(PC2,PC2save);
masse=masse{ia}(ia2,ia2);
masse = MULTIMATRIX(masse,[length(ia2),length(ia2)],[length(ia),1]);
else
%fprintf('-> Calcul de la matrice masse stochastique\n')
ind = getindices(PC);
ind2 = getindices(PC2);

PC = calc_masseuni(PC,PC2);
mpc = getmasseuni(PC);

%[repi,repj]=ind2sub([size(ind2,1),size(ind2,1)],[1:size(ind2,1)^2]');
%for k=1:getM(PC)
%CI{k}=ind2(repi,k);
%CJ{k}=ind2(repj,k);
%C{k} = [CI{k},CJ{k}];
%[ci{k},ck{k}]=find(double(mpc{k}));
%[ci{k},cj{k}]=ind2sub(size(mpc{k}),ci{k});
%ci{k} = ci{k}-1;
%cj{k} = cj{k}-1;
%ck{k} = ck{k}-1;
%c{k} = [ci{k},cj{k},ck{k}];
%end
met=4;
if met==1
    
  
masse=cell(length(PC),1);
for k=1:length(PC)
  pourcentage(k,length(PC));
  
  A = mpc{1}{ind(k,1)+1};
  masse{k} = A(ind2(:,1)+1,ind2(:,1)+1);
  for j=2:getM(PC)
  A = mpc{j}{ind(k,j)+1};
  B = A(ind2(:,j)+1,ind2(:,j)+1);
  masse{k} = masse{k} .*B;  
  end
end
%fprintf('\n')
masse = MULTIMATRIX(masse,[length(PC2),length(PC2)]);
elseif met==4
masse=cell(length(PC),1);

Kglob = [];
Lglob = [];
valglob=[];
for k=1:length(PC)
 
  pourcentage(k,length(PC),100);
  
  A = mpc{1}{ind(k,1)+1};
  B = A(ind2(:,1)+1,ind2(:,1)+1);
  [I,J,VAL] = find(B);
  K = sub2ind([length(PC2),length(PC2)],I,J);
  for j=2:getM(PC)
  A = mpc{j}{ind(k,j)+1};
  B = A(ind2(:,j)+1,ind2(:,j)+1);
  BK = B(K);
  Kt = find(BK);
  K = K(Kt);
  VAL = VAL(Kt).*BK(Kt);
  
  end
  Kglob = [Kglob;K];
  Lglob = [Lglob;repmat(k,length(K),1)];
  valglob = [valglob;VAL];
end
%fprintf('\n')
masse = sparse(Kglob,Lglob,valglob,length(PC2)*length(PC2),length(PC));
masse = MULTIMATRIX(masse,[length(PC2),length(PC2)],[length(PC),1]);
 

elseif met==2
for m=1:getM(PC)
[repig{m},repjg{m},repkg{m},valg{m}]=find(mpc{m});
end

repI = [1:length(PC)]';
repI = repmat(repI,1,length(PC));
repJ = [1:length(PC)];
repJ = repmat(repJ,length(PC),1);
repI = repI(:);
repJ = repJ(:);
indI = ind(repI,:);
indJ = ind(repJ,:);
for k=1:length(PC)
    pourcentage(k,length(PC));
  for m=1:getM(PC) 
  rep{m} = find(repkg{m}==ind(k,m)+1);
  repi{m}=repig{m}(rep{m});
  repj{m}=repjg{m}(rep{m});
  val{m}=valg{m}(rep{m});
  end
  
  
  for m=1:getM(PC)
  [temp,loc{m}]=ismember([indI(:,m),indJ(:,m)],[repi{m},repj{m}],'rows');
  
  if m==1
      f = find(temp);
  else
     f=intersect(f,find(temp));
  end
  end
  
  repIglob = repI(f);
  repJglob = repJ(f);
  m=1;
  valglob = val{m}(loc{m}(f));
  for m=2:getM(PC)
  valglob = valglob.*val{m}(loc{m}(f));    
  end
  
  
  masse{k} = sparse(repIglob,repJglob,valglob,length(PC2),length(PC2));

  
end
%fprintf('\n')

elseif met==3
for m=1:getM(PC)
[repig{m},repjg{m},repkg{m},valg{m}]=find(mpc{m});
end


for k=1:length(PC)
    pourcentage(k,length(PC));
  for m=1:getM(PC) 
  rep{m} = find(repkg{m}==ind(k,m)+1);
  repi{m}=repig{m}(rep{m});
  repj{m}=repjg{m}(rep{m});
  val{m}=valg{m}(rep{m});
  end
  repIglob = [];
  repJglob = [];
  valglob=[];
  for kk=1:length(PC2)
      m=1;
      repikk = find(repi{m}==ind(kk,m)+1);
      VALt = val{m}(repikk);
      
      [temp,temp2] = ismember(ind(:,m)+1 , repj{m}(repikk));
      repJ = find(temp);
      VAL = VALt(temp2(repJ));
      
      for m=2:getM(PC2)
         repikk = find(repi{m}==ind(kk,m)+1); 
         VALt = val{m}(repikk);
         [temp,temp2] = ismember(ind(repJ,m)+1 , repj{m}(repikk));
         repJt = find(temp);
         repJ  = repJ(repJt);
         VAL   = VAL(repJt);
         VAL = VAL.* VALt(temp2(repJt));             
      end
      
      repI = repmat(kk,[length(repJ),1]);
      repIglob = [repIglob;repI];
      repJglob = [repJglob;repJ];
      valglob = [valglob;VAL];
      
     
  end
  
  masse{k}=sparse(repIglob,repJglob,valglob,length(PC2),length(PC2));
 
end
%fprintf('\n')


end


masse = refreshsparse(masse);

PCsave = PC;
PC2save = PC2 ;


save('MULTIMASSE','masse','PCsave','PC2save');

end

PC.masse = masse;

return

tic
for p=1:M
mpc{p} = MULTIMATRIX(MPC{p});    
end
masse1=cell(PK+1,1);
for k=0:PK
  pourcentage(k,PK);
  L =[];
  value=[];
  for p=1:M
  [i,j] = find(mpc{p}{indK(k+1,p)+1});    
  l=sub2ind([P+1,P+1],i,j);
  if p>1
  L=setdiff(L,l);
  else
  L=l;
  end  
%  value = [value;temp(l)];
  end
%  masse1{k+1} = accumarray([I,J],value,[P+1,P+1],@prod); 
end
%toc
%  masse{k+1} = masse{k+1} .* MPC{j}{indK(k+1,j)+1}(ind(:,j)+1,ind(:,j)+1);    


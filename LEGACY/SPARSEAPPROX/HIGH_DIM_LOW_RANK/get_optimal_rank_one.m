function [lopt]=get_optimal_rank_one(PC,fs,rs)
Ns=size(fs,1);     
w = 1/Ns*ones(Ns,1);       
l=normalize(normalizephi(ones(1,1,PC)));
ls=cell(1,getnbgroups(PC)); 
Hs=polyval(PC,rs);
for j=1:getnbgroups(PC)
    ls{j} = Hs{j}*l{j};
end
for pp=1:30
l{0}=1;
   for i=1:getnbgroups(PC)
   l{i} =1;  
   ls{i}=ones(Ns,1);
   gammas = prod([ls{:}],2);
   l{i} = (Hs{i}'*spdiags(gammas.*w.*gammas,0,Ns,Ns)*Hs{i})\...
          (Hs{i}'*spdiags(gammas.*w,0,Ns,Ns)*fs);
   normphi = norm(l{i});
   l{i} = l{i}/normphi;
   ls{i}= Hs{i}*l{i};
   end
   l{0} = normphi;  
end
% l{1}=l{0}*l{1};
% l{0}=1;
lopt=l;
end
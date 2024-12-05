function H=purge(H,dim)
%      function A=purge(A,dim)
% Elimine les rangs = 0 en testant les dimensions dim
% Si dim non specifie, toutes les dimensions sont testees.
if nargin==1
    dim=1:totaldim(H);
end

% Reperage local :

locdim = getvar2dim(H.tree,dim,1);
unique_locdim = unique(locdim)';


for d=unique_locdim
    % Calculer le recalage :
    recalage=0;
    for dd=1:d-1
        recalage=recalage+totaldim(H.F{1,dd});
    end
    
    subdim = dim(locdim==d)-recalage;
    H.F(:,d)=cellfun(@(f) purge(f,subdim) , H.F(:,d) ,'UniformOutput',0);
    keeprank=cellfun(@(f) f.m , H.F(:,d) )~=0;
    
    H.F=H.F(keeprank,:);
    H.alpha=H.alpha(keeprank);
    H.m=length(H.alpha);
    if H.m==0
        break
    end
    
end







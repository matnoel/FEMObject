function apc=calc_ximasse(apc,PC)
% function apc=calc_ximasse(apc,PC)
% apc :  PCMATRIX
% PC : POLYCHAOS

PCa = getPC(apc);
if nargin==1
    PC = PCa;
end

if iscalcmasse(PCa)
    
    
    if (nargin == 1) || polycmp(PCa,PC)
        apc = actualise_masse(apc);
    else
        apc=actualise_masse(apc,PC);
    end
    masse=getmasse(apc);
    apc.ximasse = multimtimes(double(apc),masse);
    
    % if numel(apc)==1
    %     apc.ximasse = getmatrix(apc.ximasse);
    % end
    
else
    % warning('calcul des masse stochastique sans assembler la masse sto de POLYCHAOS')
    
    ind = getindices(PCa);
    ind2 = getindices(PC);
    
    PCa = calc_masseuni(PCa,PC);
    mpc = getmasseuni(PCa);
    
    apc.ximasse=cell(size(apc,1),size(apc,2));
    for i=1:numel(apc)
        apc.ximasse{i} = sparse(length(PC),length(PC));
    end
    
    apcvalue = double(apc);
    for k=1:length(PCa)
        
        pourcentage(k,length(PCa),100);
        A = mpc{1}{ind(k,1)+1};
        masse = A(ind2(:,1)+1,ind2(:,1)+1);
        for j=2:getM(PC)
            A = mpc{j}{ind(k,j)+1};
            B = A(ind2(:,j)+1,ind2(:,j)+1);
            masse = sparsetimes(masse,B);
        end
        
        for i=1:numel(apc)
            apc.ximasse{i} = apc.ximasse{i} +  sparsemtimes(apcvalue(i,k),masse);
        end
        
    end
    
    
end

% apc.ximasse = MULTIMATRIX(apc.ximasse);

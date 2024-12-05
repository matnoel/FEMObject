function PC = union(PC1,PC2,varargin)
% function PC = union(PC1,PC2)
% union de deux POLYCHAOS
%
% See also POLYCHAOS/intersect, POLYCHAOS/setdiff

PC.typebase = getcharin('typebase',varargin,min(PC1.typebase,PC2.typebase));

if ~polycmp(PC1.RANDPOLYS,PC2.RANDPOLYS)
    H1 = RANDPOLYS(PC1);
    H2 = RANDPOLYS(PC2);
    p1 = getorder(PC1);
    p2 = getorder(PC2);
    [ok,rep12]=ismember(H1,H2);
    rep21=find(ok);
    rep12=rep12(rep21);
    rep1 = setdiff(1:getM(H1),rep21) ;
    rep2 = setdiff(1:getM(H2),rep12) ;
    
    H12 = H1(rep21);
    if ~isempty(rep21)
        p12 = max(p1(rep21),p2(rep12));
    else
        p12=[];
    end
    PC = POLYCHAOS(RANDPOLYS(H1(rep1),H12,H2(rep2)),[p1(rep1),p12,p2(rep2)],'typebase',PC.typebase);
    
    if ischarin('independent',varargin)
        ind = getindices(PC);
        rep = sum(ind(:,1:end-1)==0,2);
        
        rep = find(rep>=getM(PC)-1);
        ind = ind(rep,:);
        
        PC.indices = ind;
        PC.P = size(PC.indices,1)-1;
    end
else
    PC = PC1;
    ind1 = getindices(PC1);
    ind2 = getindices(PC2);
    PC.indices=union(ind1(:,1:PC.M),ind2(:,1:PC.M),'rows','stable');
    switch PC.typebase
        case 1
            PC.indices(:,PC.M+1) = sum(PC.indices(:,1:PC.M),2);
        case 2
            PC.indices(:,PC.M+1) = max(PC.indices(:,1:PC.M),[],2);
    end
    if ischarin('sort',varargin)
        PC.indices = sortrows(PC.indices,size(PC.indices,2):-1:1);
    end
    PC.p = max(PC1.p,PC2.p);
    % PC.p = max(PC.indices(:,1:PC.M),[],1);
    for k=1:PC.M
        PC.n(k) = length(unique(PC.indices(:,k)));
    end
    PC.P = size(PC.indices,1)-1;
    PC.masseuni = {};
    PC.masse = {};
    PC.metricuni = {};
    PC.metric = {};
end

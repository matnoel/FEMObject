function ind=getindices(PC,varargin)
% function ind=getindices(PC)
% Get the multi-index set of POLYCHAOS PC
% 
% function ind_max=getindices(PC,'maximal')
% Get the set of maximal indices contained in the multi-index set of POLYCHAOS PC
% 
% function ind=getindices(PC,'margin')
% Get the margin of the multi-index set of POLYCHAOS PC
% 
% function ind=getindices(PC,'reducedMargin')
% Get the reduced margin of the multi-index set of POLYCHAOS PC
% 
% See also POLYCHAOS/setindices, POLYCHAOS/addindices, POLYCHAOS/removeindices

ind=PC.indices;
M=PC.M;
P=length(PC);
typebase=PC.typebase;

if ischarin('maximal',varargin)
    ind_max = [];
    for i=1:P
        neighbours = repmat(ind(i,1:M),M,1) + eye(M);
        [ok,rep]=ismember(neighbours,ind(:,1:M),'rows');
        if ~any(ok)
            ind_max = [ind_max;ind(i,1:M)];
        end
    end
    switch typebase
        case 1
            ind_max(:,M+1) = sum(ind_max(:,1:M),2);
        case 2
            ind_max(:,M+1) = max(ind_max(:,1:M),[],2);
    end
    ind = ind_max;
elseif ischarin('margin',varargin)
    ind_marg = [];
    for i=1:P
        neighbours = repmat(ind(i,1:M),M,1) + eye(M);
        ind_marg = [ind_marg;setdiff(neighbours,ind(:,1:M),'rows')];
    end
    ind_marg = unique(ind_marg,'rows');
    switch typebase
        case 1
            ind_marg(:,M+1) = sum(ind_marg(:,1:M),2);
        case 2
            ind_marg(:,M+1) = max(ind_marg(:,1:M),[],2);
    end
    ind = ind_marg;
elseif ischarin('reducedMargin',varargin)
    ind_marg = getindices(PC,'margin');
    ind_marg_red = [];
    for i=1:size(ind_marg,1)
        neighbours = repmat(ind_marg(i,1:M),M,1) - eye(M);
        [rows,cols] = find(neighbours<0);
        neighbours(rows,:) = [];
        [ok,rep]=ismember(neighbours,ind(:,1:M),'rows');
        if all(ok)
            ind_marg_red = [ind_marg_red;ind_marg(i,:)];
        end
    end
    ind = ind_marg_red;
end

if ischarin('sort',varargin)
    ind = sortrows(ind,size(ind,2):-1:1);
end

end

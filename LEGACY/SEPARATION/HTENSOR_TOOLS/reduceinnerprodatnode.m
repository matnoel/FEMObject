function M = reduceinnerprodatnode(x,y,t)
% function M = reduceinnerprodatnode(x,y,t)
%   if t is root,   M={M1,M2,[]}
%   if t is leaf,   M={[],[],M3}
%   otherwise       M={M1,M2,M3}
%   M1 corresponds to the contraction at the left children of t
%   M2 corresponds to the contraction at the right children of t
%   M3 corresponds to the contraction at the parent of t

M=cell(1,3);
p = getparents(x,t);
Mi = incomplete_innerprod(x,y,p);

%% Reduction of the tree under the node t
if ~x.is_leaf(t)
    children = x.children(t,:);
    M{1} = Mi{children(1)};
    M{2} = Mi{children(2)};
end

%% Reduction of the tree over the node t
% Contract M{ii_children} and y.B{ii} for ii in parents and ii_children
if t~=1
    Br=cell(size(Mi));
    for ii=p(2:end)
        ii_left  = x.children(ii, 1);
        ii_right = x.children(ii, 2);
        if any(ii_left==p) %
            Br{ii} = ttm(y.B{ii}, Mi{ii_right}, 2);
        else
            Br{ii} = ttm(y.B{ii}, Mi{ii_left}, 1);
        end
    end
    
    % Reduce the tree over the node t
    ps=sort(p);
    for ii=ps(1:end-1); %from top to bottom
        c=x.children(ii,:);
        ii_left=c(1);
        if any(ii_left==p)
            Q = ttt(x.B{ii},Br{ii},2,2);
            if ii>1
                if numel(Qparent) == 1
                    Q = Q*Qparent;
                else
                    if ndims(Q) == 3
                        % MATLAB automatically squeeze last dimension
                        % when its size is 1
                        Q = ttt(Q,Qparent,3,2);
                    else
                        Q = ttt(Q,Qparent,[2 4],[1 2]);
                    end
                end
            end
        else
            Q = ttt(x.B{ii},Br{ii},1,1);
            if ii>1
                if numel(Qparent) == 1
                    Q = Q*Qparent;
                else
                    if ndims(Q) == 3
                        Q = ttt(Q,Qparent,3,2);
                    else
                        Q = ttt(Q,Qparent,[2 4],[1 2]);
                    end
                end
            end
        end
        Qparent=Q;
    end
    M{3}=Q;
end

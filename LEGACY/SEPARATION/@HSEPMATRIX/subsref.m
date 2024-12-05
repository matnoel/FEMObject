function v = subsref(A,s)
% function v = subsref(A,s)

% if size(s,2)==1 && strcmp(s(1).type,'.') %% H.F / H.alpha / H.dim / H.tree / H.var...
%     switch s(1).subs
%         case 'alpha'
%             if length(s)==1
%                 v = H.alpha;
%             else
%                 v = H.alpha(s(2).subs{:});
%             end
%         case 'F'
%
%             if length(s)==1
%                 v = H.F;
%             else
%                 v = subsref(H.F,s(2:end));
%             end
%         case 'dim'
%             v = H.dim;
%         case 'm'
%             v = H.m;
%         case 'tree'
%             v = H.tree;
%         case 'var'
%             v = H.var;
%     end
%
%
% %%% Demander un alpha specifique : H(i,j)  -> cf SM
% elseif size(s,2)==1 && strcmp(s(1).type,'()')
%     if length(s.subs)==1
%         v = H.F{s.subs{1}};
%     elseif length(s.subs)==2
%         v = H.F{s.subs{1},s.subs{2}};
%     end
%
%
% end


if strcmp(s(1).type,'.')
    
    
    switch s(1).subs
        case 'alpha'
            if length(s)==1
                v = A.alpha;
            else
                v = A.alpha(s(2).subs{:});
            end
        case 'F'
            if length(s)==1
                v = A.F;
            else
                v = subsref(A.F,s(2:end));
            end
        case 'dim'
            v = A.dim;
        case 'm'
            v = A.m;
        case 'tree'
            v = A.tree;
    end
    
    
    
else
    if length(s)==1 && strcmp(s(1).type,'{}')
        if length(s.subs)==1
            v = A;
            v.F = v.F(s(1).subs{1},:);
            v.m = size(v.F);
        else
            v = A.F{s(1).subs{1},s(1).subs{2}};
        end
        
    elseif length(s)==2 && strcmp(s(1).type,'{}') && strcmp(s(2).type,'()')
        
        v = A.F{s(1).subs{1},s(1).subs{2}}(s(2).subs{:});
        
    elseif length(s)==1 && strcmp(s(1).type,'()')
        v = A.F(s(1).subs{:});
    end
end





function A = subsasgn(A,s,v)
% function A = subsasgn(A,s,v)

%
% if size(s,2)==1 && strcmp(s(1).type,'.') %% H.F / H.alpha
%
%     if strcmp(s(1).subs,'alpha')
%         if length(s)==1
%             H.alpha = v;
%         else
%             H.alpha(s(2).subs{:}) = v;
%         end
%     elseif strcmp(s(1).subs,'F')
%         if length(s)==1
%             H.F = v;
%         else
%             H = subsasgn(H,s(2:end),v);
%         end
%     elseif strcmp(s(1).subs,'m')
%         H.m=v;
%     else error('pas compris...');
%     end
%
% elseif size(s,2)>1 && strcmp(s(1).type,'{}') && strcmp(s(2).type,'()')
%
%
% end




if strcmp(s(1).type,'.')
    if strcmp(s(1).subs,'alpha')
        if length(s)==1
            A.alpha = v;
        else
            A.alpha(s(2).subs{:}) = v;
        end
    elseif strcmp(s(1).subs,'F')
        if length(s)==1
            A.F = v;
        else
            A = subsasgn(A,s(2:end),v);
        end
    end
else
    if length(s)==2
        A.F{s(1).subs{:}}(s(2).subs{:})=v;
    elseif length(s)==1
        if isa(v,'SEPMATRIX')
            A.F{s(1).subs{:}} = v;
        else
            if isa(v,'cell')
                A.F(s(1).subs{:}) = v;
            else
                A.F(s(1).subs{:}) = {v};
            end
        end
    end
    
    A.dim = size(A.F,2);
    A.m = size(A.F,1);
    
end



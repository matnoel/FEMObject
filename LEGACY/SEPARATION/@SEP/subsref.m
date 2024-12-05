function v = subsref(A,s)
% function v = subsref(A,s)

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
        otherwise
            v = eval(['A.' s(1).subs]);
    end
else
    if length(s)==1 && strcmp(s(1).type,'{}')
        if length(s.subs)==1
            if isa(s(1).subs{1},'double')
                v = A;
                v.F = v.F(s(1).subs{1},:);
                v.m = size(v.F);
            elseif isa(s(1).subs{1},'SEPMODEL')
                s(2).type = '()';
                s(2).subs = {'times'};
                v = subsref(A,s);
            end
        else
            v = A.F{s(1).subs{1},s(1).subs{2}};
        end
        
    elseif length(s)==2 && strcmp(s(1).type,'{}') && strcmp(s(2).type,'()')
        if isa(s(1).subs{1},'double')
            v = A.F{s(1).subs{1},s(1).subs{2}}(s(2).subs{:});
        elseif isa(s(1).subs{1},'SEPMODEL')
            % Operation du type S{SM}('indications')
            %             % Coder en direct
            %             if isa(s(2).subs{1},'char') && (strcmp(s(2).subs{1},'time')||strcmp(s(2).subs{1},'times'))
            %                 % u{SM}('time') : cree l'object U = masse*u
            %                 % cas SEPMATRIX{SEPMODEL}
            %                 v = A;
            %                 %SALLE :
            %                 for i=1:A.m
            %                     v.F(i,:) = cellfun(@(u,A) ...
            %                         getmatrix(switchmulti(u'*switchmulti(A))),...
            %                         v.F(i,:),getmasse(s(1).subs{1}),'uniformoutput',false);
            %                 end
            %             end
            v = SEPMATRIX(getmasse(s(1).subs{1}))*A;
        end
    elseif length(s)==1 && strcmp(s(1).type,'()')
        v = A.F(s(1).subs{:});
    end
end


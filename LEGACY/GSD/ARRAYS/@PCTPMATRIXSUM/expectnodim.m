function x = expectnodim(nodim,x,a,b)
% function x = expectnodim(nodim,x,a,b)

if nargin==2
    for i=1:length(x.funs)
        x.funs{i} = expectnodim(nodim,x.funs{i});
    end
elseif nargin==3
    if isa(x,'PCTPMATRIXSUM') && isa(a,'PCTPMATRIX')
        for i=1:length(x.funs)
            x.funs{i} = expectnodim(nodim,x.funs{i},a);
        end
    elseif isa(a,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIX')
        for i=1:length(a.funs)
            a.funs{i} = expectnodim(nodim,x,a.funs{i});
        end
        x=a;
    elseif isa(a,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIXSUM')
        z = PCTPMATRIXSUM(x.POLYCHAOSTP);
        z.funs = cell(1,length(x.funs)*length(a.funs));
        ll=0;
        for i=1:length(x.funs)
            for j=1:length(a.funs)
                ll=ll+1;
                z.funs{ll} = expectnodim(nodim,x.funs{i},a.funs{j});
            end
        end
        x=z;
    end
elseif nargin==4
    if isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIXSUM') && isa(b,'PCTPMATRIX')
        for i=1:length(a.funs)
            a.funs{i} = expectnodim(nodim,x,a.funs{i},b);
        end
        x=a;
    elseif isa(x,'PCTPMATRIXSUM') && isa(a,'PCTPMATRIX') && isa(b,'PCTPMATRIX')
        for i=1:length(x.funs)
            x.funs{i} = expectnodim(nodim,x.funs{i},a,b);
        end
    elseif isa(x,'PCTPMATRIXSUM') && isa(a,'PCTPMATRIXSUM') && isa(b,'PCTPMATRIXSUM')
        z = PCTPMATRIXSUM(x.POLYCHAOSTP);
        z.funs = cell(1,length(x.funs)*length(a.funs)*length(b.funs));
        ll=0;
        for i=1:length(x.funs)
            for j=1:length(a.funs)
                for k=1:length(b.funs)
                    ll=ll+1;
                    z.funs{ll} = expectnodim(nodim,x.funs{i},a.funs{j},b.funs{k});
                end
            end
        end
        x=z;
    elseif isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIXSUM') && isa(b,'PCTPMATRIXSUM')
        
        z = PCTPMATRIXSUM(a.POLYCHAOSTP);
        z.funs = cell(1,length(a.funs)*length(b.funs));
        ll=0;
        for j=1:length(a.funs)
            for k=1:length(b.funs)
                ll=ll+1;
                z.funs{ll} = expectnodim(nodim,x,a.funs{j},b.funs{k});
            end
        end
        x=z;
    elseif isa(x,'PCTPMATRIXSUM') && isa(a,'PCTPMATRIX') && isa(b,'PCTPMATRIXSUM')
        
        z = PCTPMATRIXSUM(x.POLYCHAOSTP);
        z.funs = cell(1,length(x.funs)*length(b.funs));
        ll=0;
        for j=1:length(x.funs)
            for k=1:length(b.funs)
                ll=ll+1;
                z.funs{ll} = expectnodim(nodim,x.funs{j},a,b.funs{k});
            end
        end
        x=z;
    else
        keyboard
        error('pas programme');
    end
else
    error('pas programme')
end


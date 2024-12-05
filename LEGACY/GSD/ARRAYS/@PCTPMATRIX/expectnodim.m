function x = expectnodim(nodim,x,a,b)
% function x = expectnodim(nodim,x,a,b)

if nargin==2
    % experance sur les dimensions autres que nodim
    rdim = setdiff(getranddim(x),nodim);
    for i=rdim
        x.phi{i} = mean(x.POLYCHAOSTP,i)'*x.phi{i};
        if isximasse(x)
            x.ximasse{i} = x.phi{i};
        end
    end
    x.isranddim(rdim)=0;
elseif nargin==3 && ~israndom(x)
    x = x*expectnodim(nodim,a);
elseif nargin==3 && ~israndom(a)
    x = expectnodim(nodim,x)*a;
elseif nargin==4 && ~israndom(x)
    x = x*expectnodim(nodim,a,b);
elseif nargin==4 && ~israndom(b)
    x = expectnodim(nodim,x,a)*b;
elseif nargin==3 && isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIX')
    x = expectnodimmtimes(nodim,x,a);
elseif nargin==4 && isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIX') && isa(b,'PCTPMATRIX')
    x = expectnodimmtimes(nodim,x,a,b);
else
    error('pas programme')
end


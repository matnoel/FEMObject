function x = expectdim(dim,x,a,b)
if ~isempty(dim)
    if nargin==2

        dim = intersect(dim,getranddim(x));
        for i=dim
            Hi = mean(x.POLYCHAOSTP,i);    
            x.phi{i} = Hi'*x.phi{i};   
        end
        if alldetdim(x)
            phi = prod(vertcat(x.phi{:}),1);
            x = x.phi0 * phi;
        end

    elseif nargin==3 && isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIX')
        if numel(a)>1 
            error('utiliser expectmtimes ou expecttimes')
        end
        dimx = intersect(dim,getranddim(x));
        dima = intersect(dim,getranddim(a));
        a = expectdim(setdiff(dima,dimx),a);
        x = expectdim(setdiff(dimx,dima),x);

        for i=dim
            x.phi{i} = a.phi{i}'*x.phi{i};     
        end
        x.phi0 = x.phi0*a.phi0;

        if alldetdim(x)
            phi = prod(vertcat(x.phi{:}),1);
            x = x.phi0 * phi;
        end


    elseif nargin==4 && isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIX') && isa(b,'PCTPMATRIX')

        if numel(a)>1 || numel(b)>1
            error('utiliser expectmtimes ou expecttimes')
        end
        dimx = intersect(dim,getranddim(x));
        dima = intersect(dim,getranddim(a));
        dimb = intersect(dim,getranddim(a));
        a = expectdim(setdiff(dima,union(dimx,dimb)),a);
        b = expectdim(setdiff(dimb,union(dimx,dima)),b);
        x = expectdim(setdiff(dimx,union(dima,dimb)),x);
        dimaxb = unique([dima,dimb,dimx]);
        x.phi0 = a.phi0 * x.phi0 * b.phi0 ; 
        for i=dimaxb
            ximasse = get_ximasse(x,i);   
            x.phi{i} = a.phi{i}'*ximasse*b.phi{i};   
        end
        for i=setdiff(dimaxb,dimx)
            ximasse = get_ximasse(x,i);   
            x.phi{i} = a.phi{i}'*ximasse*b.phi{i};      
        end

        phi = prod(vertcat(x.phi{:}),1);
        x = (a.phi0*x.phi0*b.phi0)* phi;

    else
        error('pas programme')
    end
end

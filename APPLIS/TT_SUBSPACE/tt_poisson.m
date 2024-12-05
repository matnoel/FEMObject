classdef tt_poisson
    properties
        A
        b
        Att
        btt
    end
    
    
    methods
        function obj = tt_poisson(d,n)
            D = DOMAIN(1);
            S = cell(d,1);
            Kb = BILINFORM(1,1);
            Mb = BILINFORM(0,0);
            As = cell(d,1);
            for mu = 1:d
                S{mu} = mesh(D,n(mu));
                S{mu} = createddlnode(S{mu},DDL('u'));
                S{mu} = addcl(S{mu},POINT([0;1]),'u');
                for i = 1:d
                    if mu == i
                        As{mu}{i} = Kb{S{mu}}(:,:);
                    else
                        As{mu}{i} = Mb{S{mu}}(:,:);
                    end
                end
            end
            As = TSPACE_OPERATORS(As);
            Aa = CANONICAL_CORE(ones(d,1),d);
            obj.A = LRTENSOR(Aa,As);
            obj.Att = LRTENSOR(TT_CORE(Aa),As);
            
            f = LINFORM(0,1);
            bs = cell(d,1);
            for mu = 1:d
                bs{mu} = full(f{S{mu}}(:));
            end
            bs = TSPACE_VECTORS(bs);
            bc = CANONICAL_CORE(1,d);
            obj.b = LRTENSOR(bc,bs);
            obj.btt = LRTENSOR(TT_CORE(bc),bs);
        end
    end
    
    
end

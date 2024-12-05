function c = simplify(c,varargin)

repsimp = zeros(1,length(c.tensorfuns));
onedim = zeros(1,length(c.tensorfuns));
V = PRODUCTSPACE();
    for i=1:length(c.tensorfuns)
        c.tensorfuns{i} = simplify(c.tensorfuns{i});
        if isa(c.tensorfuns{i},'double')
        repsimp(i) = 1;
        elseif isa(c.tensorfuns{i},'TENSORPRODUCTSPACE')          
        V = union(V,getproducspace(c.tensorfuns{i}));
        end
    end
c.PRODUCTSPACE = V;

if all(repsimp)
   z = 0;
   for i=1:length(c.tensorfuns)
   z = z + c.tensorfuns{i};
   end
   c=z;
end
    
function Mfe=findpolyfe(H)

Mfe = [];
for k=1:length(H.h)
if isa(H.h{k},'POLYFE')
 Mfe = [Mfe,k];
end
end

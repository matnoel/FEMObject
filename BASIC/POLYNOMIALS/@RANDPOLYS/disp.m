function disp(polys)
% function disp(polys)

for k=1:polys.M
    if ~isempty(polys.h{k})
        fprintf('{%d} ',k)
        disp(polys.h{k})
    end
end

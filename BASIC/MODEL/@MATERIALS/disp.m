function disp(mats)
% function disp(mats)

for k=1:mats.n
    if ~isempty(mats.MAT{k})
        fprintf('{%d} ',k)
        disp(mats.MAT{k})
    end
end

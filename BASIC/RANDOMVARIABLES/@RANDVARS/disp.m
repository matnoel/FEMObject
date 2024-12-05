function disp(rvs)
% function disp(rvs)

for k=1:rvs.M
    if ~isempty(rvs.RV{k})
        fprintf('{%d} ',k)
        disp(rvs.RV{k})
    end
end

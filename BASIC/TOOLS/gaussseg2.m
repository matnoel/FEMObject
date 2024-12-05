function gauss = gaussseg2(order)
% function gauss = gaussseg2(order)

n = ceil((order+1)/2);
try
    load('gaussseg2');
    if length(gaussorder)>=n && ~isempty(gaussorder{n})
        gauss = gaussorder{n};
        return
    end
    if length(gaussorder)<n
        gaussorder = [gaussorder , cell(1,n-length(gaussorder))] ;
    end
catch
    gaussorder = cell(1,n); 
end

gauss = calc_gausspoints(POLYLEGENDRE(),n);
gauss.w = gauss.w*2;
gaussorder{n} = gauss;

save('gaussseg2','gaussorder');

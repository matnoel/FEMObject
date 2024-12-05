function x = transfer(apc,x)
% function x = transfer(PC,x)
% Transfer matrix sample x on random variables associated to PC

if isa(apc,'PCMODEL') || isa(apc,'FESTOMODEL') || isa(apc,'PCMATRIX') || isa(apc,'PCARRAY')
    % warning(['The function to be decomposed is a function of random variables associated to ' class(apc)])
    if isa(apc,'PCMODEL') || isa(apc,'FESTOMODEL')
        apc = vertcat(apc{:});
    elseif isa(apc,'PCMATRIX')
        apc = apc(:);
    end
    x = randomeval(apc,x);
    x = double(x)';
end

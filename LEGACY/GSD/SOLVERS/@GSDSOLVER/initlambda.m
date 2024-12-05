function l0 = initlambda(GSD,PC)
% function l0 = initlambda(GSD,PC)

switch getparam(GSD,'inittype')
    case 'one'
        l0 = one(PC);
    case 'random'
        l0 = rand(PC);
    case 'allone'
        l0 = ones(length(PC),1);
        l0 = PCMATRIX(l0,[1,1],PC);
end

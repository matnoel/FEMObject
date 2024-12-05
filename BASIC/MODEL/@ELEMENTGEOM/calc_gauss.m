function gauss = calc_gauss(elem,choice,varargin)
% function gauss = calc_gauss(elem,choice,varargin)

if getparam(elem,'initializeBN')
    gauss = getparam(elem,'gauss');
    return
end

p = calc_gaussorder(elem,choice);

gauss = elem_gauss(elem,p);
gauss = permutegaussND(gauss);

%gauss.coord = permute(gauss.coord,[4,2,3,1]);
%gauss.w     = permute(gauss.w(:),[2,3,4,1]);



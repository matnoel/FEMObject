function source = randomizeSource(model,maxSourceValue)
% source = randomizeSource(model,maxSourceValue)
% Generate random source. Based on formatSource.m.

switch randi(3)
    case 1 % double
        nbCellDoF = getNbCellDoF(model) ;
        nbDomainDoF = getNbDomainDoF(model) ;
        sourceSize = {1 ;... % uniform source
            [nbCellDoF 1] ;... % scalar field over cell
            [nbDomainDoF 1]} ; % field over domain
        source = maxSourceValue*rand(sourceSize{randi(numel(sourceSize))}) ;
    case 2 % function handle
        source = @(x) maxSourceValue*rand(length(x),1) ;
    case 3 % char
        source = {'corrector1';'corrector2'} ;
        source = source{randi(numel(source))} ;
end

end
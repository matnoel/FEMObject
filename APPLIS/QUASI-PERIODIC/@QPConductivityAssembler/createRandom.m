function Assembler = createRandom(varargin)
% Assembler = createRandom(varargin)

% Hard-coded bounds
maxFieldNb = 9 ;
maxFieldValue = 100 ;
minFieldValue = 1 ;


% Parse varargin and and random values to missing properties
if ~ischarin('model',varargin)
    varargin = setcharin('model',varargin,...
        QPModel.createRandom(varargin{:})) ;
end
if ~ischarin('fields',varargin) && ~ischarin('patterns',varargin)
    fieldNb = randi(maxFieldNb) ;
    model = getcharin('model',varargin) ;
    fields = maxFieldValue*rand(getNbCellDoF(model),fieldNb) ;
    fields = max(fields,minFieldValue) ;
    varargin = setcharin('fields',varargin,fields) ;
else
    fieldNb = size(getcharin('fields',varargin),2) ;
end
if ~ischarin('distributor',varargin)
    if fieldNb == 1
        distributor = @(n) {(1:n)'} ;
    else
        proba = getcharin('probability',varargin,...
            rand(1,fieldNb-1)/(fieldNb-1)) ;
        if numel(proba)<fieldNb
            proba = [1-sum(proba) proba] ;
        end
        distributor = @(n) dealMultinomial(proba,n) ;
    end
    varargin = setcharin('distributor',varargin,distributor) ;
end

Assembler = QPConductivityAssembler(varargin{:}) ;
end
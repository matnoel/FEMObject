function [feModel,mesherTime] = buildFEModel(pb,method)
% [feModel,mesherTime] = buildFEModel(pb,method)

if nargin == 1
    method = 1 ;
end

[feModel,mesherTime] = buildFEModel(pb.model,method) ;

end
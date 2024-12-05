function [rep,reps] = isparamin(P,varargin)
% function rep = isparamin(P,'paramname1','paramname2')
% verifient que les parametres  'paramname1','paramname2' sont bien des paramtres de P

reps = zeros(1,length(varargin)) ;
for i=1:length(varargin)
    if isfield(P.param,varargin{i})
        reps(i) = 1;
    end
end

rep = all(reps);

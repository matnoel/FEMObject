function varargout = checkparamname(P,varargin)
% function rep = checkparamname(P,'paramname1','paramname2')
% verifient que les parametres  'paramname1','paramname2' sont bien des paramtres de P

[rep,reps] = isparamin(P,varargin{:});

if ~rep
    fprintf('Les parametres suivants ne sont pas dans l''object PARAMETERS\n')
    disp(varargin(find(reps==0)));
    error(' ')
end


if nargout==1
    varargout{1} = rep;
end


if nargout==2
    varargout{2} = reps;
end


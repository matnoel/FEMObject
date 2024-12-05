function [b,varargout]=pcsystemupdate(b,varargin)
noupdate = ischarin('noupdate',varargin);
if noupdate
varargin = delcharin('noupdate',varargin);
end
if israndom(b)
PC = getPC(b);
else
    for i=1:length(varargin)
    if israndom(varargin{i})
        PC = getPC(varargin{i});
    break
    end
    end
if ~exist('PC')
    error('le systeme n''est pas aleatoire')
end
end


varargout=cell(1,length(varargin)+1);

for i=1:length(varargin)
varargout{i} = varargin{i};
if ~noupdate
if isa(varargout{i},'TIMEMATRIX') && israndom(varargout{i}) 
varargout{i}.value = calc_masse(varargout{i}.value,PC);
elseif israndom(varargout{i})
varargout{i} = calc_masse(varargout{i},PC);    
end
end
end

varargout{end}=PC;
function  plot(T,f,varargin)
% function plot(T,f,varargin)
% T : TIMEMODEL
% f : 'double' representant une ou plusieurs fonctions temporelles
% varargin : arguments pour la fonction plot
% appel de plot(t,f,varargin{:})
% t pas de temps
%

if isa(f,'TIMEMODEL') && ~isa(f,'TIMEMATRIX')
    plot(f,T,varargin{:})
else
    
    if ~isa(f,'TIMEMATRIX')
        f = TIMEMATRIX(f,T);
    end
    
    plot(f,varargin{:});
end
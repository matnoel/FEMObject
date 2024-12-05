function [u,result] = dsolve(GSD,T,varargin)
% function u = dsolve(GSD,T,b,A,B,varargin)


if ~isa(T,'TIMEMODEL') || strcmp(class(T),'TIMEMODEL')
    error('rentrer un SOLVER TEMPOREL')
end

if israndom(varargin)
    
switch getparam(GSD,'type')    
    case 'power'
     [u,result] = dsolve_power_gsd_sto(GSD,T,varargin{:})   ;
    case 'arnoldi'
     [u,result] = dsolve_arnoldi_gsd_sto(GSD,T,varargin{:})   ;  
end

else

switch getparam(GSD,'type')    
    case 'power'
     [u,result] = dsolve_power_gsd_time(GSD,T,varargin{:})   ;
    case 'arnoldi'
     [u,result] = dsolve_arnoldi_gsd_time(GSD,T,varargin{:})   ;  
end

    
end
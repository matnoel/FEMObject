function [u,result] = newtonglob(GSD,varargin)
% function [u,result] = newtonglob(GSD,N,S,b,X,varargin)
% GSD : GSDSOLVER
% N : NEWTONSOLVER
% S : MODEL
% b : second membre (PCMATRIX) ...
% X : POLYCHAOS

switch getparam(GSD,'type')
    case 'power'
        [u,result] = newtonglob_power(GSD,varargin{:})   ;
    case 'arnoldi'
        [u,result] = newtonglob_arnoldi(GSD,varargin{:})   ;
    otherwise
        error('pas prevu')
        
end

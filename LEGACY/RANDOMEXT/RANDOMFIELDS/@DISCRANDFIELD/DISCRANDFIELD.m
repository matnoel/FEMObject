function DRF = DISCRANDFIELD(varargin)
% function RF = DISCRANDFIELD(RANDFIELD,PCRADIALMATRIX,MODEL)
% MODEL
% PCRADIALMATRIX with object of type FENODEFIELD
% RANDFIELD
% function RF = DISCRANDFIELD(RANDFIELD,PCARRAY,MODEL)

DRF.RF = getclassin('RANDFIELD',varargin);
PCR = getclassin('PCRADIALMATRIX',varargin);
if isempty(PCR)
    PCR = getclassin('PCMATRIX',varargin);
end
if isempty(PCR)
    error('donner un PCRADIALMATRIX ou PCMATRIX') 
end
DRF.PCR = PCR;    
DRF.S = getclassin('MODEL',varargin);

DRF = class(DRF,'DISCRANDFIELD');

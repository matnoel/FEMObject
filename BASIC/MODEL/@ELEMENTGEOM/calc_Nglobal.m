function [N,detJ,x]=calc_Nglobal(elem,xnode,xgauss,varargin)

[N,detJ,x]=calc_N(elem,xnode,xgauss,varargin{:});


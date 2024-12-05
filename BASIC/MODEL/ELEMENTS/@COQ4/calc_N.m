function [N,varargout] = calc_N(elem,xnode,xgauss)
% function [N,varargout] = calc_N(elem,xnode,xgauss)

varargout = cell(1,nargout-1);
[N,varargout{:}] = calc_N(elem.QUA4,xnode,xgauss);
P = calc_P(elem);
N = N*P;

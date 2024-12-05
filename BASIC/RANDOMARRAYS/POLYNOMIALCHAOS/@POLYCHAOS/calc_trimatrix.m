function t=calc_trimatrix(PC,a,varargin)
% function t=calc_trimatrix(PC,a)
% A completer si mapping spatial, donne par a.k

try
    t=getmasse(PC);
catch
    t=getmasse(calc_masse(PC));
end

function A = calc_matrix(a,S,varargin)
% function A = calc_matrix(a,S,varargin)

if getn(a)~=2
    error('il faut une BILINFORM')
else
    A = calc_matrix(a.forms{1},S,varargin{:});
    for i=2:length(a.forms)
        A = A + calc_matrix(a.forms{i},S,varargin{:});
    end
end


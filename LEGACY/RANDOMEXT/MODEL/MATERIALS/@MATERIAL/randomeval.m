function mat = randomeval(mat,varargin)
% function mat = randomeval(mat,varargin)

 mat.param = funrandomparam(mat.param,@randomeval,varargin{:});

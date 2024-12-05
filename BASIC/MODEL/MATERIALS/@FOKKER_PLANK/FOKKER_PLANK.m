function mat = FOKKER_PLANK(varargin)
% function mat = FOKKER_PLANK('D1',D1,'D2',D2)
% d(D1 u) - d2(D2 u)

mat = struct();

param = struct(varargin{:});

matp = MATERIAL('FOKKER_PLANK',param);
mat = class(mat,'FOKKER_PLANK',matp);

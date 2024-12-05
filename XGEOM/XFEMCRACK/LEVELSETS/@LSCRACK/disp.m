function disp(c,varargin)
% function disp(c,varargin)

fprintf(' CRACK SUPPORT ')
disp(c.LEVELSETS{1},varargin{:})
for k=1:getnbtip(c)
    fprintf(' CRACK TIP #%d',k)
    disp(c.LEVELSETS{1+k},varargin{:})
end





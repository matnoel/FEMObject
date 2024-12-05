function disp(PC,varargin)
% function disp(PC,varargin)

v = struct('M',PC.M,'p',PC.p,'n',PC.n,'P',PC.P,'typebase',PC.typebase);

display(struct(v))
if ~ischarin('nopolys',varargin)
    disp(PC.RANDPOLYS)
end

% if ~isempty(PC.germ)
%     fprintf('Germ')
%     disp(PC.germ)
% end
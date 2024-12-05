function hix = subsref(h,s)
% function hix = subsref(h,s)

switch length(s.subs)
    case 2
        hix = polyval(h,s.subs{1},s.subs{2});
    case 1
        hix = polycoeff(h,s.subs{1});
    case 0
        error('no subsref')
end